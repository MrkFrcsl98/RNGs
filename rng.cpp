#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <algorithm>
#include <array>
#include <cerrno>
#include <cstring>
#include <chrono>
#include <thread>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>


namespace RNG
{
/**
 * @brief Advanced Pseudo-Random Number Generator Implementation
 *
 * This class implements four distinct PRNG algorithms, each serving different purposes:
 *
 * 1. Linear Congruential Generator (LCG):
 *    - Fastest but statistically weakest
 *    - Uses linear recurrence: X(n+1) = (a * X(n) + c) mod m
 *    - Good for simple applications where speed is critical
 *
 * 2. Xorshift Generator (XSG):
 *    - Uses bitwise operations for speed
 *    - Better statistical properties than LCG
 *    - Good balance of speed and quality
 *
 * 3. Mersenne Twister (MT19937):
 *    - Industry standard for statistical applications
 *    - Extremely long period of 2^19937-1
 *    - Passes most statistical tests
 *
 * 4. Permuted Congruential Generator (PCG):
 *    - Modern algorithm combining LCG with permutation
 *    - Excellent statistical properties
 *    - Small state size but high-quality output
 */
class PRNG
{
  private:
    /**
     * Mersenne Twister Constants:
     *
     * N = 624: State size, determines the period (2^19937-1)
     * M = 397: Middle word offset for twist transformation
     * MATRIX_A: Twist transformation constant (0x9908b0df)
     *          - Used in the twist operation to ensure maximum period
     *          - Represents a primitive polynomial over GF(2)
     *
     * UPPER_MASK: Most significant w-r bits (0x80000000)
     *            - Used to extract the high bits during twist
     *            - Represents the most significant bit in a 32-bit word
     *
     * LOWER_MASK: Least significant r bits (0x7fffffff)
     *            - Used to extract the low bits during twist
     *            - Represents all bits except the most significant
     */
    static constexpr std::size_t N = 624;
    static constexpr std::size_t M = 397;
    static constexpr std::size_t MATRIX_A = 0x9908b0dfUL;
    static constexpr std::size_t UPPER_MASK = 0x80000000UL;
    static constexpr std::size_t LOWER_MASK = 0x7fffffffUL;

    /**
     * LCG and PCG Constants:
     *
     * MULTIPLIER (6364136223846793005ULL):
     * - Carefully chosen prime number
     * - Ensures maximum period in the LCG sequence
     * - Properties:
     *   1. Large enough to provide good distribution
     *   2. Congruent to 5 modulo 8 for optimal lattice structure
     *   3. Passes spectral test for dimensions up to 8
     *
     * INCREMENT (1442695040888963407ULL):
     * - Must be odd for full period
     * - Helps avoid correlations between successive values
     * - Properties:
     *   1. Large prime number
     *   2. No common factors with modulus
     *
     * MODULUS (2^32):
     * - Power of 2 for efficient computation
     * - Allows use of bit operations instead of division
     * - Provides period of 2^32 for LCG
     */
    static constexpr std::size_t MULTIPLIER = 6364136223846793005ULL;
    static constexpr std::size_t INCREMENT = 1442695040888963407ULL;
    static constexpr std::size_t MODULUS = (1ULL << 32);

    /**
     * State Variables:
     *
     * state: Current internal state of the generator
     *        - Updated with each number generation
     *        - Must never be zero for LCG
     *
     * inc: PCG stream identifier
     *      - Must be odd
     *      - Different values produce different sequences
     *
     * mt: Mersenne Twister state array
     *     - Holds 624 32-bit values
     *     - Updated in blocks for efficiency
     *
     * mti: Current position in Mersenne Twister array
     *      - Ranges from 0 to N-1
     *      - Triggers state update when reaches N
     */
    std::size_t state;
    std::size_t inc;
    std::array<std::size_t, N> mt;
    int mti;

    /**
     * @brief Initialize Mersenne Twister state array using advanced recurrence formula
     *
     * Detailed initialization process:
     * 1. Set first element directly from seed
     * 2. Each subsequent element is generated using the formula:
     *    mt[i] = f * (mt[i-1] xor (mt[i-1] >> 30)) + i
     *
     * The formula components:
     * - f = 1812433253: Multiplier constant (empirically derived)
     * - Right shift by 30: Removes lower order bits
     * - XOR operation: Mixes bits for better randomization
     * - Addition of index i: Ensures unique values across array
     *
     * This initialization ensures:
     * - High quality initial state
     * - No zero elements in the state array
     * - Good bit distribution from the start
     *
     * @param seed Initial random seed value
     */
    void init_mersenne_twister(std::size_t seed)
    {
        mt[0] = seed;
        for (mti = 1; mti < N; mti++)
        {
            // Complex initialization formula broken down:
            // 1. Take previous state: mt[mti-1]
            // 2. Right shift by 30 bits: mt[mti-1] >> 30
            // 3. XOR with original: mt[mti-1] ^ (mt[mti-1] >> 30)
            // 4. Multiply by constant: 1812433253UL * (...)
            // 5. Add index for uniqueness: + mti
            mt[mti] = (1812433253UL * (mt[mti - 1] ^ (mt[mti - 1] >> 30)) + mti);
        }
    }

  public:
    /**
     * @brief Initialize all PRNG components with careful state setup
     *
     * Constructor operations:
     * 1. Sets initial state from seed
     * 2. Creates odd increment for PCG (sequence << 1 | 1)
     * 3. Initializes MT array using complex recurrence
     *
     * Parameter details:
     * @param seed: Initial entropy source
     *             - Defaults to current time
     *             - Should be unpredictable for security
     *
     * @param sequence: PCG stream identifier
     *                 - Must be unique for different streams
     *                 - Converted to odd number internally
     */
    PRNG(std::size_t seed = std::time(nullptr), std::size_t sequence = 1)
        : state(seed), inc((sequence << 1) | 1) // Ensure odd increment
          ,
          mti(N)
    {
        init_mersenne_twister(seed);
    }

    /**
     * @brief Generate a random number using Linear Congruential Generator
     *
     * Formula: state = (MULTIPLIER * state + INCREMENT) % MODULUS
     *
     * @param min Minimum value (inclusive)
     * @param max Maximum value (inclusive)
     * @return std::size_t Random number in range [min, max]
     * @throws std::invalid_argument if min >= max
     */
    std::size_t LCG(const std::size_t min, const std::size_t max)
    {
        // Check if the provided minimum value is greater than or equal to the maximum value
        // If so, throw an exception to indicate that the arguments are invalid.
        if (min >= max)
            throw std::invalid_argument("min must be less than max");

        // Update the state using a linear congruential generator formula.
        // MULTIPLIER, INCREMENT, and MODULUS are constants that define the behavior of the LCG.
        // This formula generates a new state based on the previous state.
        state = (MULTIPLIER * state + INCREMENT) % MODULUS;

        // Calculate a random number in the range [min, max].
        // The expression (state % (max - min + 1)) gives a value in the range [0, max - min],
        // and adding min shifts it to the desired range [min, max].
        return min + (state % (max - min + 1));
    }

    /**
     * @brief Generate a random number using Xorshift Generator
     *
     * Uses a series of shift and XOR operations to generate random numbers
     *
     * @param min Minimum value (inclusive)
     * @param max Maximum value (inclusive)
     * @return std::size_t Random number in range [min, max]
     * @throws std::invalid_argument if min >= max
     */
    std::size_t XSG(const std::size_t min, const std::size_t max)
    {
        // Check if the provided minimum value is greater than or equal to the maximum value
        // If so, throw an exception to indicate that the arguments are invalid.
        if (min >= max)
            throw std::invalid_argument("min must be less than max");

        // Perform a series of bitwise operations on the 'state' variable to generate new
        // pseudorandom bits. These operations are part of a PRNG algorithm, which mixes the bits of
        // 'state' to produce a new value that is more uniformly distributed.
        state ^= state << 13; // XOR the state with itself shifted left by 13 bits
        state ^= state >> 7;  // XOR the state with itself shifted right by 7 bits
        state ^= state << 17; // XOR the state with itself shifted left by 17 bits

        // Calculate a random number in the range [min, max].
        // The expression (state % (max - min + 1)) gives a value in the range [0, max - min],
        // and adding min shifts it to the desired range [min, max].
        return min + (state % (max - min + 1));
        // Check if the provided minimum value is greater than or equal to the maximum value
        // If so, throw an exception to indicate that the arguments are invalid.
        if (min >= max)
            throw std::invalid_argument("min must be less than max");

        // Perform a series of bitwise operations on the 'state' variable to generate new
        // pseudorandom bits. These operations are part of a PRNG algorithm, which mixes the bits of
        // 'state' to produce a new value that is more uniformly distributed.
        state ^= state << 13; // XOR the state with itself shifted left by 13 bits
        state ^= state >> 7;  // XOR the state with itself shifted right by 7 bits
        state ^= state << 17; // XOR the state with itself shifted left by 17 bits

        // Calculate a random number in the range [min, max].
        // The expression (state % (max - min + 1)) gives a value in the range [0, max - min],
        // and adding min shifts it to the desired range [min, max].
        return min + (state % (max - min + 1));
    }
    /**
     * @brief Mersenne Twister implementation with detailed bit manipulation
     *
     * Algorithm steps:
     * 1. State Array Generation (when mti >= N):
     *    - Process array in blocks for efficiency
     *    - Apply twist transformation:
     *      a. Combine upper bits of current word with lower bits of next
     *      b. Right shift combined value
     *      c. XOR with magic constant if lowest bit is 1
     *      d. XOR result with word M positions ahead
     *
     * 2. Tempering Transform (improves statistical properties):
     *    y ^= (y >> 11)                    // Right shift tempering
     *    y ^= (y << 7) & 0x9d2c5680UL     // Left shift tempering with mask
     *    y ^= (y << 15) & 0xefc60000UL    // Additional left shift tempering
     *    y ^= (y >> 18)                    // Final right shift tempering
     *
     * Each tempering step:
     * - Improves randomness of output
     * - Breaks up patterns in raw state
     * - Uses carefully chosen constants
     *
     * @param min Lower bound (inclusive)
     * @param max Upper bound (inclusive)
     * @return Tempered random number in [min, max]
     */
    std::size_t MersenneTwister(const std::size_t min, const std::size_t max)
    {
        // Check if the provided minimum value is greater than or equal to the maximum value.
        // If this condition is true, throw an invalid argument exception to prevent further
        // processing.
        if (min >= max)
            throw std::invalid_argument("min must be less than max");

        // Declare a variable to hold the generated random number.
        std::size_t y;

        // Define a static array for the Mersenne Twister algorithm.
        // This array is used to help with the generation of new random numbers.
        // The first element is 0, and the second element is the constant MATRIX_A.
        static const std::size_t mag01[2] = {0x0UL, MATRIX_A};

        // Check if the index for the Mersenne Twister state array (mti) has reached the size of the
        // array (N). If mti is greater than or equal to N, it means we need to generate a new batch
        // of random numbers.
        if (mti >= N)
        {
            int kk; // Declare a loop variable for iterating through the state array.

            // Generate N words at a time using the Mersenne Twister algorithm.
            // This loop fills the state array with new values based on the previous values.
            for (kk = 0; kk < N - M; kk++)
            {
                // Combine the upper and lower bits of the current and next state values.
                y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);

                // Update the current state value using the Mersenne Twister recurrence relation.
                // The new value is derived from the value M positions ahead and the combined bits.
                mt[kk] = mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x1UL];
            }

            // Continue filling the state array for the remaining elements.
            for (; kk < N - 1; kk++)
            {
                // Again combine the upper and lower bits of the current and next state values.
                y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);

                // Update the current state value using the recurrence relation.
                mt[kk] = mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
            }

            // Handle the last element of the state array separately to complete the generation.
            y = (mt[N - 1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
            mt[N - 1] = mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

            // Reset the index to 0, indicating that we have filled the state array.
            mti = 0;
        }

        // Retrieve the next random number from the state array using the current index (mti).
        y = mt[mti++];

        // Apply the tempering process to the generated number.
        // Tempering is a technique used to improve the statistical properties of the generated
        // numbers. It modifies the bits of the number to ensure better randomness.
        y ^= (y >> 11);                // Shift right by 11 bits and XOR with the original number.
        y ^= (y << 7) & 0x9d2c5680UL;  // Shift left by 7 bits, AND with a mask, and XOR.
        y ^= (y << 15) & 0xefc60000UL; // Shift left by 15 bits, AND with a mask, and XOR.
        y ^= (y >> 18);                // Shift right by 18 bits and XOR with the original number.

        // Return a random number within the specified range [min, max].
        // The result is calculated by taking the modulus of y with the range size and adding min.
        return min + (y % (max - min + 1));
    }

    /**
     * @brief PCG algorithm core with detailed permutation steps
     *
     * Operation sequence:
     * 1. LCG State Update:
     *    state = state * MULTIPLIER + inc
     *    - Multiplier and increment ensure full period
     *    - State change is reversible
     *
     * 2. Output Permutation:
     *    a. Right shift state by 18 bits
     *    b. XOR with original state
     *    c. Right shift result by 27 bits
     *    d. Calculate rotation amount (state >> 59)
     *    e. Perform rotate right operation
     *
     * Permutation benefits:
     * - Improves statistical properties
     * - Makes output non-linear
     * - Harder to predict than raw LCG
     *
     * @return Permuted random value
     */
    std::size_t PCGenerator()
    {
        // Store the current state in a temporary variable 'oldstate'
        std::size_t oldstate = state;

        // Update the state using a linear congruential generator formula
        // MULTIPLIER is a constant that determines how the state evolves,
        // and 'inc' is an increment value that can be used to introduce variability.
        state = oldstate * MULTIPLIER + inc;

        // Perform a Xorshift operation on the old state
        // This operation shifts the bits of 'oldstate' and applies XOR to create a new value.
        std::size_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;

        // Calculate the rotation amount by shifting 'oldstate' to the right by 59 bits
        // This determines how many bits to rotate the 'xorshifted' value.
        std::size_t rot = oldstate >> 59u;

        // Return the final random value by rotating the 'xorshifted' value
        // The rotation is done using bitwise operations. The result combines
        // the shifted bits to produce a uniformly distributed random number.
        return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    }

    /**
     * @brief Generate a random number using PCG
     *
     * Uses the PCG algorithm with improved statistical properties
     *
     * @param min Minimum value (inclusive)
     * @param max Maximum value (inclusive)
     * @return std::size_t Random number in range [min, max]
     * @throws std::invalid_argument if min >= max
     */
    std::size_t generatePC(const std::size_t min, const std::size_t max)
    {
        // Check if the provided minimum value is greater than or equal to the maximum value
        // If so, throw an exception to indicate that the arguments are invalid.
        if (min >= max)
            throw std::invalid_argument("min must be less than max");

        // Calculate the range of possible values, which is the difference between max and min,
        // plus one to include both endpoints.
        std::size_t range = max - min + 1;

        // Calculate the threshold value to ensure uniform distribution.
        // This threshold is used to discard values that would lead to bias in the random number
        // generation.
        std::size_t threshold = -range % range;

        // Infinite loop to repeatedly generate random numbers until a valid one is found.
        while (true)
        {
            // Generate a random number using a pseudorandom number generator (PCGenerator).
            std::size_t r = PCGenerator();

            // Check if the generated random number is greater than or equal to the threshold.
            // If it is, we can safely map it to the desired range without bias.
            if (r >= threshold)
            {
                // Return a random number in the range [min, max].
                // The expression (r % range) gives a value in the range [0, range - 1],
                // and adding min shifts it to the desired range.
                return min + (r % range);
            }
        }
    }

    /**
     * @brief Reseed all generators with a new seed value
     *
     * @param new_seed New seed value
     */
    void reseed(std::size_t new_seed)
    {
        state = new_seed;
        init_mersenne_twister(new_seed);
    }
};

// macro defined for CSPRNG output byte size
#ifndef __DEFAULT_OUTPUT_BYTE_SIZE__
#define __DEFAULT_OUTPUT_BYTE_SIZE__ 32
#endif

class CSPRNG
{
  public:
    explicit CSPRNG() noexcept {};

    ~CSPRNG() noexcept {};

    const std::size_t generate(const std::size_t min, const std::size_t max)
    {
        int bytes = 0;
        try
        {
            this->_entropy_collect_source(); // collect entropy
            if (this->_sEntropy.urandom.empty())
                return 0;
            std::size_t state{0}; // initialize state
            for (auto c : this->_sEntropy.urandom)
            {
                state += (int)c; // add hex byte as int to the state
            }
            state += state + std::time(nullptr) / 2; // add more entropy to state

            // apply XORshift
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            bytes = min + (state % (max - min + 1)); // store into bytes
        }
        catch (const std::exception &e)
        {
            std::cerr << "CSPRNG Error: " << e.what() << "\n";
        }
        return bytes;
    };

  private:
    std::string _hex_bytes;
    struct entropySource
    {
        std::string urandom{};
        std::string cpu{};   // reserve for future
        std::string mem{};   // reserved for future
        std::string procs{}; // reserved for future
    } _sEntropy;

    void _entropy_collect_source()
    {
        this->_read_urandom();
    };

    void _read_urandom() __attribute__((hot, stack_protect))
    {
        std::ifstream FILE("/dev/urandom", std::ios::binary);
        if (!FILE.is_open())
        {
            std::string errmsg("Error Code: ");
            errmsg += strerror(errno);
            errmsg += ", Cannot Open File /dev/urandom";
            throw std::runtime_error("Cannot open file /dev/urandom !!");
        }
        const unsigned short int THRESHOLD = __DEFAULT_OUTPUT_BYTE_SIZE__; // read 32 bytes
        unsigned short int byte_count = 0; // keep track of bytes read
        std::array<unsigned char, THRESHOLD> fileContent;
        std::string hexBytes;
        char bytes_read; // current char to read in recursion
        // read from file, not more than threshold
        FILE.read(reinterpret_cast<char *>(fileContent.data()), THRESHOLD);
        bytes_read = FILE.gcount();
        hexBytes.reserve(fileContent.size() * 2); // reserve space for converter content
        if (fileContent.size() > 0 && fileContent.size() <= THRESHOLD)
        {
            for (unsigned short int i = 0; i < bytes_read; ++i)
            {
                this->_sEntropy.urandom += fileContent[i];
            }
        }
    };
};
}; // namespace RNG


int main(int argc, char** argv) {

    // Not Cryptographically Secure PRNGs
    RNG::PRNG prng;
    std::cout << "Linear Congruential Generator:   " << prng.LCG(10, 10000) << "\n";
    std::cout << "XORShift Generator:              " << prng.XSG(10, 10000) << "\n";
    std::cout << "Permuted Congruential Generator: " << prng.generatePC(10, 10000) << "\n";
    std::cout << "Mersenne Twister Generator:      " << prng.MersenneTwister(10, 10000) << "\n";

    // Cryptographically Secure PRNG
    RNG::CSPRNG csprng;
    std::cout << "CryptoSecPRNG Generator:         " << csprng.generate(10, 10000) << "\n";

    

    return 0;
}
