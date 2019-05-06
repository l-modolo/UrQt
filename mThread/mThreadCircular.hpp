/*
UrQt quality and poly nucleotide trimming tool
Copyright (C) 2013  Laurent Modolo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DEF_mThreadCircular
#define DEF_mThreadCircular

#include <exception>
#include <stdexcept>
#include <iostream>
#include <mutex>

using namespace std;

template <typename T>
class mThreadCircular
{
  public:
  mThreadCircular();
  ~mThreadCircular();
  void init(int size);

  void add(T* x, int number);
  T* pop();
  T* pop(int* number);
  inline bool can_add(int number);
  inline bool can_get();
  inline int size();
  void print();

  private:
  bool C_init;
  int C_size;
  T** C_loop;
  int* C_number;
  int C_filled;
  int C_start;
  int C_stop;
  int C_expected_start_number;
  mutex C_onebyone;
};

// initialization of the Circular container
template<typename T>
mThreadCircular<T>::mThreadCircular()
{
  try
  {
    C_init = false;
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : mThreadCircular<T>::mThreadCircular()" << endl;
  }
}

// initialization of the Circular container
template<typename T>
mThreadCircular<T>::~mThreadCircular()
{
  try
  {
    delete[] C_loop;
    delete[] C_number;
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : mThreadCircular<T>::~mThreadCircular()" << endl;
  }
}

template<typename T>
void mThreadCircular<T>::init(int size)
{
  try
  {
    C_filled = 0;
    C_size = size;
    C_loop = nullptr;
    C_number = nullptr;
    C_loop = new T*[C_size];
    C_number = new int[C_size];
    for(int i = 0; i < C_size; i++)
    {
      C_loop[i] = nullptr;
      C_number[i] = -1;
    }
    C_start = 0;
    C_stop = 0;
    C_expected_start_number = 0;
    C_init = true;
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : void mThreadCircular<T>::init(int size)" << endl;
  }
}

template<typename T>
void mThreadCircular<T>::add(T* x, int number)
{
  try
  {
    unique_lock<mutex> lk(C_onebyone);
    if(!C_init)
      throw logic_error("mThreadCircular not initisialized");

    int pos = number%C_size;
    if(C_number[pos] != -1)
    {
      throw logic_error("mThreadCircular full");
    }
    C_loop[pos] = x;
    C_number[pos] = number;
    if(number > C_number[C_stop])
      C_stop = pos;
    C_filled++;
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : void mThreadCircular<T>::add(T* x, int number)" << endl;
    exit(-1);
  }
}

template<typename T>
T* mThreadCircular<T>::pop()
{
  try
  {
    unique_lock<mutex> lk(C_onebyone);
    if(!C_init)
      throw logic_error("mThreadCircular not initisialized");

    int pos = C_start;
    if(pos == -1)
      return nullptr;

    C_start++;
    if(C_size <= C_start)
      C_start = C_start - C_size;
    C_expected_start_number = C_number[pos]+1;
    T* task = C_loop[pos];

    C_number[pos] = -1;
    C_loop[pos] = nullptr;
    C_filled--;
    return task;
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : T* mThreadCircular<T>::pop()" << endl;
    exit(-1);
  }
}

template<typename T>
T* mThreadCircular<T>::pop(int* number)
{
  try
  {
    unique_lock<mutex> lk(C_onebyone);
    if(!C_init)
      throw logic_error("mThreadCircular not initisialized");

    int pos = C_start;
    if(pos == -1)
      return nullptr;

    C_start++;
    if(C_size <= C_start)
      C_start = C_start - C_size;

    C_expected_start_number = C_number[pos]+1;
    T* task = C_loop[pos];
    *number  = C_number[pos];

    C_number[pos] = -1;
    C_loop[pos] = nullptr;
    C_filled--;
    return task;
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : T* mThreadCircular<T>::pop()" << endl;
    exit(-1);
  }
}

template<typename T>
inline bool mThreadCircular<T>::can_add(int number)
{
  try
  {
    unique_lock<mutex> lk(C_onebyone);
    if(!C_init)
      throw logic_error("mThreadCircular not initisialized");

    int pos = number%C_size;
    if(C_number[pos] == -1 && number - C_expected_start_number < C_size)
      return true;
    else
      return false;
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : bool mThreadCircular<T>::can_add()" << endl;
  }
}

template<typename T>
inline bool mThreadCircular<T>::can_get()
{
  try
  {
    unique_lock<mutex> lk(C_onebyone);
    if(!C_init)
      throw logic_error("mThreadCircular not initisialized");

    if(C_filled > 0 && C_number[C_start] != -1)
      return true;
    return false;
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : bool mThreadCircular<T>::can_get()" << endl;
  }
}

template<typename T>
inline int mThreadCircular<T>::size()
{
  try
  {
    unique_lock<mutex> lk(C_onebyone);
    if(!C_init)
      throw logic_error("mThreadCircular not initisialized");
    return C_filled;
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : int mThreadCircular<T>::size()" << endl;
  }
}


template<typename T>
void mThreadCircular<T>::print()
{
  try
  {
    unique_lock<mutex> lk(C_onebyone);
    if(!C_init)
      throw logic_error("mThreadCircular not initisialized");
    for(int i = 0; i < C_size; i++)
    {
      if(i == C_start)
        cout << "[";
      cout << C_number[i];
      if(i == C_stop)
        cout << "]";
      cout << " ";
    }
    cout << " | " << C_expected_start_number << endl;
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : int mThreadCircular<T>::size()" << endl;
  }
}



#endif


