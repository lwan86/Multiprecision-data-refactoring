#ifndef _MDR_REFACTOR_INTERFACE_HPP
#define _MDR_REFACTOR_INTERFACE_HPP

namespace MDR {
    namespace concepts {

        // refactor: a general interface for scnetific data refactor
        template<class T>
        class RefactorInterface {
        public:

            virtual ~RefactorInterface() = default;

            virtual void refactor(T const * data_, const std::vector<uint32_t>& dims, uint8_t target_level, uint8_t num_bitplanes) = 0;

            virtual uint8_t * write_metadata(uint32_t& size) const = 0;

            virtual uint8_t * get_data(const std::vector<double>& eb, std::vector<int>& positions) = 0;

            virtual void print() const = 0;
        };
    }
}
#endif
