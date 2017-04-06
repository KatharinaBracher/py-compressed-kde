#pragma once

#include <vector>
#include <algorithm>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

static const double SQRT2 = std::sqrt(2);

typedef double value; // double or float

// Fast approximations of power, exponential and logarithm
#define cast_uint32_t (uint32_t)

// these function require float types!!
float fastpow2 (float p);
float fastexp (float p);
float fastlog2 (float x);
float fastlog (float x);


value circular_difference(value a, value b);

size_t hash_combine( size_t lhs, size_t rhs );

void multiply_add_vectors( const std::vector<std::vector<value>> & v, unsigned int N, value weight, value * result );
void multiply_add_vectors( const std::vector<std::vector<value>> & v, unsigned int N, value weight, value * result, const std::vector<bool> & valid );
void add_assign_vectors( const std::vector<std::vector<value>> & v, unsigned int N, value weight, value * result );
void add_assign_vectors( const std::vector<std::vector<value>> & v, unsigned int N, value weight, value * result, const std::vector<bool> & valid );


template <typename T>
bool isunique( std::vector<T> v ) {
    
    // vector is passed by value, so we can manipulate it here
    std::sort( v.begin(), v.end() );
    return (std::adjacent_find( v.begin(), v.end() ) == v.end());
    
}

