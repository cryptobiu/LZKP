//
// Created by roee on 9/13/18.
//

#ifndef __LZKP_UTILS_H_FILE__
#define __LZKP_UTILS_H_FILE__


namespace lzkp {


template <typename T>
T** allocate2D(unsigned int N, unsigned int M) {
  T** ptr = new T*[N]; // Allocate pointers
  T* pool = new T[N * M]; // Allocate pool

  for (unsigned i = 0; i < N; ++i, pool += M)
    ptr[i] = pool;

  return ptr;
}

template <typename T>
void free2Darray(T **ptr) {
  if (ptr) {
    delete[] ptr[0]; // Free the pool
    delete[] ptr;
  }
}

template <typename T>
T* allocate1D(unsigned int N) {
  T* ptr = new T[N]; // Allocate pointers

  return ptr;
}

template <typename T>
void free1Darray(T *ptr) {
  if (ptr) {
    delete[] ptr;
  }
}


}


#endif // __LZKP_UTILS_H_FILE__