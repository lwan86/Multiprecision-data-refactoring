#ifndef _MDR_PERBIT_BP_ENCODER_HPP
#define _MDR_PERBIT_BP_ENCODER_HPP

#include "BitplaneEncoderInterface.hpp"
#include <bitset>
namespace MDR {
    class BitEncoder{
    public:
        BitEncoder(uint64_t * stream_begin_pos){
            stream_begin = stream_begin_pos;
            stream_pos = stream_begin;
            buffer = 0;
            position = 0;
        }
        void encode(uint64_t b){
            buffer += b << position;
            position ++;
            if(position == 64){
                *(stream_pos ++) = buffer;
                buffer = 0;
                position = 0;
            }
        }
        void flush(){
            if(position){
                *(stream_pos ++) = buffer;
                buffer = 0;
                position = 0;
            }
        }
        uint32_t size(){
            return (stream_pos - stream_begin);
        }
    private:
        uint64_t buffer = 0;
        uint8_t position = 0;
        uint64_t * stream_pos = NULL;
        uint64_t * stream_begin = NULL;
    };

    class BitDecoder{
    public:
        BitDecoder(uint64_t const * stream_begin_pos){
            stream_begin = stream_begin_pos;
            stream_pos = stream_begin;
            buffer = *(stream_pos ++);
            position = 64;
        }
        uint8_t decode(){
            uint8_t b = buffer & 1u;
            buffer >>= 1;
            position --;
            if(position == 0){
                buffer = *(stream_pos ++);
                position = 64;
            }
            return b;
        }
        uint32_t size(){
            return (stream_pos - stream_begin);
        }
    private:
        uint64_t buffer = 0;
        uint8_t position = 0;
        uint64_t const * stream_pos = NULL;
        uint64_t const * stream_begin = NULL;
    };

    #define PER_BIT_BLOCK_SIZE 1
    // per bit bitplane encoder that encodes data by bit using T_stream type buffer
    template<class T_data, class T_stream>
    class PerBitBPEncoder : public concepts::BitplaneEncoderInterface<T_data> {
    public:
        PerBitBPEncoder(){
            static_assert(std::is_floating_point<T_data>::value, "PerBitBPEncoder: input data must be floating points.");
            static_assert(!std::is_same<T_data, long double>::value, "PerBitBPEncoder: long double is not supported.");
            static_assert(std::is_unsigned<T_stream>::value, "PerBitBPEncoder: streams must be unsigned integers.");
            static_assert(std::is_integral<T_stream>::value, "PerBitBPEncoder: streams must be unsigned integers.");
        }

        std::vector<uint8_t *> encode(T_data const * data, int32_t n, int32_t exp, uint8_t num_bitplanes, std::vector<uint32_t>& stream_sizes) const {
            assert(num_bitplanes > 0);
            // determine block size based on bitplane integer type
            const int32_t block_size = PER_BIT_BLOCK_SIZE;
            std::vector<uint8_t> starting_bitplanes = std::vector<uint8_t>((n - 1)/block_size + 1, 0);
            stream_sizes = std::vector<uint32_t>(num_bitplanes, 0);
            // define fixed point type
            using T_fp = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;
            std::vector<uint8_t *> streams;
            for(int i=0; i<num_bitplanes; i++){
                streams.push_back((uint8_t *) malloc(2 * n / UINT8_BITS + sizeof(uint64_t)));
            }
            std::vector<BitEncoder> encoders;
            for(int i=0; i<streams.size(); i++){
                encoders.push_back(BitEncoder(reinterpret_cast<uint64_t*>(streams[i])));
            }
            T_data const * data_pos = data;
            for(int i=0; i<n - block_size; i+=block_size){
                T_stream sign_bitplane = 0;
                for(int j=0; j<block_size; j++){
                    T_data cur_data = *(data_pos++);
                    T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
                    bool sign = cur_data < 0;
                    int64_t fix_point = (int64_t) shifted_data;
                    T_fp fp_data = sign ? -fix_point : +fix_point;
                    bool first_bit = true;
                    for(int k=num_bitplanes - 1; k>=0; k--){
                        uint8_t index = num_bitplanes - 1 - k;
                        uint8_t bit = (fp_data >> k) & 1u;
                        encoders[index].encode(bit);
                        if(bit && first_bit){
                            encoders[index].encode(sign);
                            first_bit = false;
                        }
                    }                    
                }
            }
            // leftover
            {
                int rest_size = n % block_size;
                if(rest_size == 0) rest_size = block_size;
                for(int j=0; j<rest_size; j++){
                    T_data cur_data = *(data_pos++);
                    T_data shifted_data = ldexp(cur_data, num_bitplanes - exp);
                    bool sign = cur_data < 0;
                    int64_t fix_point = (int64_t) shifted_data;
                    T_fp fp_data = sign ? -fix_point : +fix_point;
                    bool first_bit = true;
                    for(int k=num_bitplanes - 1; k>=0; k--){
                        uint8_t index = num_bitplanes - 1 - k;
                        uint8_t bit = (fp_data >> k) & 1u;
                        encoders[index].encode(bit);
                        if(bit && first_bit){
                            encoders[index].encode(sign);
                            first_bit = false;
                        }
                    }                    
                }
            }
            for(int i=0; i<num_bitplanes; i++){
                encoders[i].flush();
                stream_sizes[i] = encoders[i].size() * sizeof(uint64_t);
            }
            return streams;
        }

        T_data * decode(const std::vector<uint8_t const *>& streams, int32_t n, int exp, uint8_t num_bitplanes) const {
            assert(num_bitplanes > 0);
            const int32_t block_size = PER_BIT_BLOCK_SIZE;
            // define fixed point type
            using T_fp = typename std::conditional<std::is_same<T_data, double>::value, uint64_t, uint32_t>::type;
            T_data * data = (T_data *) malloc(n * sizeof(T_data));
            std::vector<BitDecoder> decoders;
            for(int i=0; i<streams.size(); i++){
                decoders.push_back(BitDecoder(reinterpret_cast<uint64_t const*>(streams[i])));
                decoders[i].size();
            }
            // decode
            T_data * data_pos = data;
            for(int i=0; i<n - block_size; i+=block_size){
                for(int j=0; j<block_size; j++){
                    T_fp fp_data = 0;
                    // decode each bit of the data for each level component
                    bool first_bit = true;
                    bool sign = false;
                    for(int k=num_bitplanes - 1; k>=0; k--){
                        uint8_t index = num_bitplanes - 1 - k;
                        uint8_t bit = decoders[index].decode();
                        fp_data += bit << k;
                        if(bit && first_bit){
                            // decode sign
                            sign = decoders[index].decode();
                            first_bit = false;
                        }
                    }
                    T_data cur_data = ldexp((T_data)fp_data, - num_bitplanes + exp);
                    *(data_pos++) = sign ? -cur_data : cur_data;
                }
            }
            // leftover
            {
                int rest_size = n % block_size;
                if(rest_size == 0) rest_size = block_size;
                for(int j=0; j<rest_size; j++){
                    T_fp fp_data = 0;
                    // decode each bit of the data for each level component
                    bool first_bit = true;
                    bool sign = false;
                    for(int k=num_bitplanes - 1; k>=0; k--){
                        uint8_t index = num_bitplanes - 1 - k;
                        uint8_t bit = decoders[index].decode();
                        fp_data += bit << k;
                        if(bit && first_bit){
                            // decode sign
                            sign = decoders[index].decode();
                            first_bit = false;
                        }
                    }
                    T_data cur_data = ldexp((T_data)fp_data, - num_bitplanes + exp);
                    *(data_pos++) = sign ? -cur_data : cur_data;
                }
            }
            return data;
        }

        void print() const {
            std::cout << "Per-bit bitplane encoder" << std::endl;
        }

    };
}
#endif
