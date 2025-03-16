# RNGs

# Overview
- PRNG
- - Linear Congruential Generator
  - XORShift Generator
  - Permuted Congruential Generator
  - Mersenne Twister Generator

- CSPRNG
- - Generator
 
# Usage 

### Linear Congruential Generator(LCG)
```cpp
RNG::PRNG prng;
std::uint64_t gen = prng.LCG(min, max);
```

### Permuted Congruential Generator(PCG)
```cpp
RNG::PRNG prng;
std::uint64_t gen = prng.generatePC(min, max);
```

### XOR Shift Generator(XSG)
```cpp
RNG::PRNG prng;
std::uint64_t gen = prng.XSG(min, max);
```

### Mersenne Twister Generator(MTG)
```cpp
RNG::PRNG prng;
std::uint64_t gen = prng.MersenneTwister(min, max);
```
