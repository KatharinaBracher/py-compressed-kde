#include "common.hpp"
#include <stdint.h>
#include <stdexcept>

float fastpow2 (float p) {
  float offset = (p < 0) ? 1.0f : 0.0f;
  float clipp = (p < -126) ? -126.0f : p;
  int w = clipp;
  float z = clipp - w + offset;
  union { uint32_t i; float f; } v = { cast_uint32_t ( (1 << 23) * (clipp + 121.2740575f + 27.7280233f / (4.84252568f - z) - 1.49012907f * z) ) };

  return v.f;
}

float fastexp (float p) {
  return fastpow2 (1.442695040f * p);
}

float fastlog2 (float x) {
  union { float f; uint32_t i; } vx = { x };
  union { uint32_t i; float f; } mx = { (vx.i & 0x007FFFFF) | 0x3f000000 };
  float y = vx.i;
  y *= 1.1920928955078125e-7f;

  return y - 124.22551499f
           - 1.498030302f * mx.f 
           - 1.72587999f / (0.3520887068f + mx.f);
}

float fastlog (float x) {
  return 0.69314718f * fastlog2 (x);
}

value circular_difference(value a, value b) {
    value tmp = (b-a);
    if (tmp<0) {tmp=-tmp;}
    tmp = M_PI - tmp;
    if (tmp<0) {tmp=-tmp;}
    tmp = M_PI - tmp;
    return tmp;
}


size_t hash_combine( size_t lhs, size_t rhs ) {
    lhs^= rhs + 0x9e3779b9 + (lhs<<6) + (lhs>>2);
    return lhs;
}
    
    
void multiply_add_vectors( const std::vector<std::vector<value>> & v, unsigned int N, value weight, value * result ) {
    
    unsigned int D = v.size();
    
    std::vector< std::vector<value>::const_iterator > initer(D);
    std::vector< std::vector<value>::const_iterator > inend(D);
    for (unsigned int d=0;d<D;++d) {
        initer[d] = v[d].cbegin();
        inend[d] = v[d].cend();
    }
    
    for (unsigned int k=0; k<N; ++k) {
        *result++ += std::accumulate( initer.cbegin(), initer.cend(), weight, []( const value & a, const std::vector<value>::const_iterator & b) { return a * *b;} );
               
        for (int d=D-1; d>=0; --d) {
            ++initer[d];
            if (initer[d]>=inend[d]) {
                initer[d] = v[d].cbegin();
            } else {
                break;
            }
        }
        
    }
    
}


void multiply_add_vectors( const std::vector<std::vector<value>> & v, unsigned int N, value weight, value * result, const std::vector<bool> & valid ) {
    
    unsigned int D = v.size();
    
    if (valid.size()!=N) {
        throw std::runtime_error("Internal error: validity vector size and grid size do not match.");
    }
    
    std::vector< std::vector<value>::const_iterator > initer(D);
    std::vector< std::vector<value>::const_iterator > inend(D);
    for (unsigned int d=0;d<D;++d) {
        initer[d] = v[d].cbegin();
        inend[d] = v[d].cend();
    }
    
    for (unsigned int k=0; k<N; ++k) {
        
        if (valid[k]) {
            *result += std::accumulate( initer.cbegin(), initer.cend(), weight, []( const value & a, const std::vector<value>::const_iterator & b) { return a * *b;} );
        }
        
        ++result;
        
        for (int d=D-1; d>=0; --d) {
            ++initer[d];
            if (initer[d]>=inend[d]) {
                initer[d] = v[d].cbegin();
            } else {
                break;
            }
        }
        
    }
    
}


void add_assign_vectors( const std::vector<std::vector<value>> & v, unsigned int N, value factor, value * result ) {
    
    unsigned int D = v.size();
    
    std::vector< std::vector<value>::const_iterator > initer(D);
    std::vector< std::vector<value>::const_iterator > inend(D);
    for (unsigned int d=0;d<D;++d) {
        initer[d] = v[d].cbegin();
        inend[d] = v[d].cend();
    }
    
    for (unsigned int k=0; k<N; ++k) {
        *result++ = std::accumulate( initer.cbegin(), initer.cend(), factor, []( const value & a, const std::vector<value>::const_iterator & b) { return a + *b;} );
               
        for (int d=D-1; d>=0; --d) {
            ++initer[d];
            if (initer[d]>=inend[d]) {
                initer[d] = v[d].cbegin();
            } else {
                break;
            }
        }
        
    }
    
}


void add_assign_vectors( const std::vector<std::vector<value>> & v, unsigned int N, value factor, value * result, const std::vector<bool> & valid ) {
    
    unsigned int D = v.size();
    
    if (valid.size()!=N) {
        throw std::runtime_error("Internal error: validity vector size and grid size do not match.");
    }
    
    std::vector< std::vector<value>::const_iterator > initer(D);
    std::vector< std::vector<value>::const_iterator > inend(D);
    for (unsigned int d=0;d<D;++d) {
        initer[d] = v[d].cbegin();
        inend[d] = v[d].cend();
    }
    
    for (unsigned int k=0; k<N; ++k) {
        
        if (valid[k]) {
            *result = std::accumulate( initer.cbegin(), initer.cend(), factor, []( const value & a, const std::vector<value>::const_iterator & b) { return a + *b;} );
        }
        
        ++result;
        
        for (int d=D-1; d>=0; --d) {
            ++initer[d];
            if (initer[d]>=inend[d]) {
                initer[d] = v[d].cbegin();
            } else {
                break;
            }
        }
        
    }
    
}
