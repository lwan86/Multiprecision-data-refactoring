#include <iostream>
#include <cstdlib>
#include <vector>
#include <dirent.h>

#include <adios2.h>

#include "../include/Decomposer/Decomposer.hpp"


bool has_suffix(const std::string &str, const std::string &suffix)
{
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

int main(int argc, char *argv[])
{
    std::string inputDirName;
    std::string outputFileName;
    std::string outputEngine;
    std::string varType;
    size_t varDimensions;
    std::vector<size_t> varShape;

    for (size_t i = 0; i < argc; i++)
    {
        std::string arg = argv[i];
        if (arg == "-i" || arg == "--input")
        {
            if (i+1 < argc)
            {
                inputDirName = argv[i+1];
            }
            else
            {
                std::cerr << "--input option requires one argument." << std::endl;
                return 1;
            }            
        }
        else if (arg == "-t" || arg == "--type")
        {
            if (i+1 < argc)
            {
                varType = argv[i+1];
            }
            else
            {
                std::cerr << "--type option requires one argument." << std::endl;
                return 1;
            }            
        }
        else if (arg == "-d" || arg == "--dimension")
        {
            if (i+1 < argc)
            {
                varDimensions = atoi(argv[i+1]);
            }
            else
            {
                std::cerr << "--dimension option requires [# of dimensions] to be set first." << std::endl;
                return 1;
            }  
            if (varDimensions)
            {
                if (i+1+varDimensions < argc)
                {
                    for (size_t j = i+2; j < i+2+varDimensions; j++)
                    {
                        varShape.push_back(atoi(argv[j]));
                    }
                    
                }
                
            }             
        }
        else if (arg == "-o" || arg == "output")
        {
            if (i+2 < argc)
            {
                outputFileName = argv[i+1];
                outputEngine = argv[i+2]; 
            }
            else
            {
                std::cerr << "--output option requires two arguments." << std::endl;
                return 1;
            }            
        } 
    }

    for (size_t i = 0; i < varShape.size(); i++)
    {
        std::cout << varShape[i] << " ";
    }
    std::cout << std::endl;

    adios2::ADIOS adios;
    adios2::IO writer_io = adios.DeclareIO("WriterIO");
    adios2::Engine writer_engine = writer_io.Open(outputFileName, adios2::Mode::Write); 
    
    if (auto dir = opendir(inputDirName.c_str())) {
        while (auto f = readdir(dir)) {
            if (f->d_name[0] == '.')
                continue; // Skip everything that starts with a dot

            //printf("File: %s\n", f->d_name);
            std::string fileName(f->d_name);
            //std::cout << fileName << std::endl;
            if (has_suffix(fileName, ".dat"))
            {
                //std::cout << fileName << std::endl;
                if (!inputDirName.empty() && inputDirName.back() != '/')
                {
                    inputDirName += '/';
                }
                std::string fullPath = inputDirName+fileName;
                std::string varName = fileName.substr(0, fileName.find("."));
                std::cout << varName << std::endl;
                if (varType == "float")
                {
                    using T = float;
                    size_t n_elements = 0;
                    auto inputData = MGARD::readfile<T>(fullPath.c_str(), n_elements);
                    std::cout << n_elements << std::endl;
                    std::vector<size_t> varStart(varDimensions, 0); 
                    std::vector<size_t> varCount(varShape);
                    adios2::Variable<T> rawVariable = writer_io.DefineVariable<T>(varName, varShape, varStart, varCount);
                    writer_engine.Put(rawVariable, inputData.data(), adios2::Mode::Sync);
                }
                
            }
            
        }
        closedir(dir);
    }    
    writer_engine.Close();
}

