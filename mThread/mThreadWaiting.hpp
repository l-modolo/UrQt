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

#ifndef DEF_mThreadWaiting
#define DEF_mThreadWaiting

#include <future>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <string>
#include "mThreadCircular.hpp"

using namespace std;

template <typename T>
class mThreadWaiting
{
	public:
	mThreadWaiting(int size);
	~mThreadWaiting();
	
	void add(T* x);
	void get(T** x, int* number);
	void get(T*** x, int* number, int job_to_get);
	
	void stop();

	string output();
	
	protected:
	
	inline void push_back(T* x);
	inline T* pop_front(int* number);
	
	inline int can_add(int number);
	inline int can_get();
	// int size();
	inline int true_size();
	int mThreadWaiting_iteration;
	
	// communication between the thread
	mutex mThreadWaiting_onebyone;
	mutex mThreadWaiting_empty; // if the list is empty we wait new jobs
	mutex mThreadWaiting_full; // if the list is full we wait before new jobs
	condition_variable mThreadWaiting_empty_cond;
	condition_variable mThreadWaiting_full_cond;
	
	// status of the waiting list
	int mThreadWaiting_size; // size of the list
	int mThreadWaiting_pos_front; // we isrun jobs at the front
	int mThreadWaiting_pos_back; // we add new jobs at the back
	// deque<T*> mThreadWaiting_waiting; // mThreadWaiting_iteration of thread waiting
	mThreadCircular<T> mThreadWaiting_loop;

	bool mThreadWaiting_isrun;
};

template <typename T>
void mThreadWaiting<T>::stop()
{
	unique_lock<mutex> lk(mThreadWaiting_onebyone);
	mThreadWaiting_isrun = false;
	mThreadWaiting_empty_cond.notify_one();
}

// initialisation of the wainting list
template<typename T>
mThreadWaiting<T>::mThreadWaiting(int size)
{
	// try
	// {
		mThreadWaiting_isrun = true;
		mThreadWaiting_size = size;
		mThreadWaiting_pos_front = -1;
		mThreadWaiting_pos_back = -1;
		// mThreadWaiting_waiting.clear();
		mThreadWaiting_loop.init(size);
		mThreadWaiting_iteration = -1;
	// }
	// catch(exception const& e)
	// {
	// 	cerr << "ERROR : " << e.what() << " in : mThreadWaiting<T>::mThreadWaiting(int size)" << endl;
	// }
}

template<typename T>
mThreadWaiting<T>::~mThreadWaiting()
{
	unique_lock<mutex> full(mThreadWaiting_full);
	while(true_size() >= 1)
		mThreadWaiting_empty_cond.wait(full);
}

// we add jobs at the waiting list
template<typename T>
void mThreadWaiting<T>::add(T* x)
{
	// try
	// {
		unique_lock<mutex> full(mThreadWaiting_full);
		while(!can_add(mThreadWaiting_iteration+1))
			mThreadWaiting_full_cond.wait(full); // if the list is full we wait
		push_back(x); // when there is room we add the job
	// }
	// catch(exception const& e)
	// {
	// 	cerr << "ERROR : " << e.what() << " in : void mThreadWaiting<T>::add(T const & x)" << endl;
	// 	exit(-1);
	// }
}

// we load the next jobs
template<typename T>
void mThreadWaiting<T>::get(T** x, int* number)
{
	// try
	// {
		unique_lock<mutex> empty(mThreadWaiting_empty);
		while(!can_get())
			mThreadWaiting_empty_cond.wait(empty); // if the list in empty we wait
		*x = pop_front(number); // when there is a job we get it
	// }
	// catch(exception const& e)
	// {
	// 	cerr << "ERROR : " << e.what() << " in : T mThreadWaiting<T>::get()" << endl;
	// 	exit(-1);
	// }
}

// we load the job_to_get next jobs
template<typename T>
void mThreadWaiting<T>::get(T*** x, int* number, int job_to_get)
{
	// try
	// {
		unique_lock<mutex> empty(mThreadWaiting_empty);
		for(int i; i < job_to_get; i++)
		{
			while(!can_get())
				mThreadWaiting_empty_cond.wait(empty); // if the list in empty we wait
			*x = pop_front(number); // when there is a job we get it
		}
	// }
	// catch(exception const& e)
	// {
	// 	cerr << "ERROR : " << e.what() << " in : T mThreadWaiting<T>::get()" << endl;
	// 	exit(-1);
	// }
}

// adding a job to the list
template<typename T>
inline void mThreadWaiting<T>::push_back(T* x)
{
	unique_lock<mutex> lk(mThreadWaiting_onebyone); // we add one job at the time
	mThreadWaiting_iteration++;
	mThreadWaiting_loop.add(x, mThreadWaiting_iteration);
	mThreadWaiting_empty_cond.notify_one(); // we signal one thread waiting for a job that the list is not empty anymore (in the get() function)
}

// getting a job from the list
template<typename T>
inline T* mThreadWaiting<T>::pop_front(int* number)
{
	// try
	// {
		unique_lock<mutex> lk(mThreadWaiting_onebyone); // we get one job at the time
		T* value = nullptr;

		value = mThreadWaiting_loop.pop(number);

		if(!mThreadWaiting_isrun) // if we get the stop signal there will be no new job added so their is no reason for the thread to wait for the list to fill up
			mThreadWaiting_empty_cond.notify_all();
		else // else we allow only one thread to continue
			mThreadWaiting_full_cond.notify_one();
		return value;
	// }
	// catch(exception const& e)
	// {
	// 	cerr << "ERROR : " << e.what() << " in : T mThreadWaiting<T>::pop_front()" << endl;
	// 	exit(-1);
	// }
}

template<typename T>
inline int mThreadWaiting<T>::can_add(int number)
{
	// try
	// {
		unique_lock<mutex> lk(mThreadWaiting_onebyone);
		if(mThreadWaiting_isrun) // if we didn't get the stop signal we proceed
			return mThreadWaiting_loop.can_add(number);
		else // if we got the stop signal we can finish every thing even if the list is full
			throw logic_error("try to add task after the stop signal");
	// }
	// catch(exception const& e)
	// {
	// 	cerr << "ERROR : " << e.what() << " in : int mThreadWaiting<T>::can_add(int number)" << endl;
	// 	exit(-1);
	// }
}

template<typename T>
inline int mThreadWaiting<T>::can_get()
{
	// try
	// {
		unique_lock<mutex> lk(mThreadWaiting_onebyone);
		if(mThreadWaiting_isrun)
			return mThreadWaiting_loop.can_get();
		return true;
	// }
	// catch(exception const& e)
	// {
	// 	cerr << "ERROR : " << e.what() << " in : int mThreadWaiting<T>::can_get()" << endl;
	// 	exit(-1);
	// }
}

template<typename T>
inline int mThreadWaiting<T>::true_size()
{
	unique_lock<mutex> lk(mThreadWaiting_onebyone);
	return mThreadWaiting_loop.size();
}

template<typename T>
string mThreadWaiting<T>::output()
{
		return "W : "+to_string(mThreadWaiting_loop.size());
}

#endif

