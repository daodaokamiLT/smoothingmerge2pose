#pragma once

#include <vector>
#include <atomic>
#include <mutex>
namespace svv_fusion{

    // dataset size is 0~capacity, has capacity+1 size but save capacity elements.
    template <class T>
    class CircleQue{
        private:
            int capacity_;
            std::mutex dmtx_;
            std::vector<T> datas_;
            std::atomic_int start_;
            std::atomic_int end_;
            T tmp_ = T();
        public:
            CircleQue(const int sz = 100):capacity_(sz){
                datas_.resize(capacity_+1);
                start_.store(0);
                end_.store(0);
            }
            int size(){
                if(empty())
                    return 0;
                if(full())
                    return capacity_;
                if(end_ == start_){
                    printf("error, empty() judgment error....\n");
                    exit(-1);
                }
                else if(end_ > start_)
                    return end_ - start_;
                else
                    return end_ + capacity_+1 - start_;
            }
            
            int capacity(){
                return capacity_;
            }
            bool empty(){
                if(start_ == end_){
                    return true;
                }
                return false;
            }

            bool full(){
                if(start_ == 0 && end_ == capacity_)
                    return true;
                if(start_ - 1 == end_ || end_ + 1 == start_)
                    return true;
                return false;
            }

            bool push_front(T& data){
                dmtx_.lock();
                if(full()){
                    return false;
                }
                else{
                    --start_;
                    if(start_ == -1){
                        start_ = capacity_;
                    }
                    datas_[start_] = data;
                }
                dmtx_.unlock();
            }

            void push_front_focus(T& data){
                if(full()){
                    dmtx_.lock();
                    --end_;
                    --start_;
                    if(end_ == -1){
                        end_ = capacity_;
                    }
                    if(start_ == -1){
                        start_ = capacity_;
                    }
                    datas_[start_] = data;
                    dmtx_.unlock();
                }
                else{
                    push_front(data);
                }
            }

            bool push_back(T& data){
                if(full()){
                    return false;
                }
                else{
                    dmtx_.lock();
                    datas_[end_] = data;
                    ++end_;
                    if(end_ == capacity_+1){
                        end_ = 0;
                    }
                    dmtx_.unlock();
                }
                return true;
            }

            void push_back_focus(T& data){
                if(full()){
                    dmtx_.lock();
                    datas_[end_] = data;
                    ++end_;
                    ++start_;
                    if(end_ == capacity_+1){
                        end_ = 0;
                    }
                    if(start_ == capacity_+1){
                        start_ = 0;
                    }
                    dmtx_.unlock();
                }
                else{
                    push_back(data);
                }
            }

            bool pop_front(){
                if(empty()){
                    return false;
                }

                dmtx_.lock();
                ++start_;
                if (start_ == capacity_+1)
                {
                    start_ = 0;
                }
                dmtx_.unlock();
                return true;
            }

            bool pop_back(){

                if(empty()){
                    return false;
                }
                dmtx_.lock();
                --end_;
                if(end_ == -1){
                    end_ = capacity_;
                }
                dmtx_.unlock();
            }

            bool index(size_t i, T& data){
                if(i >= 0 || i < size()){
                    dmtx_.lock();
                    int index = 0;
                    index = start_+i;
                    if(start_+i >= capacity_+1){
                        index = start_ + i - capacity_ - 1;
                    }
                    data = datas_[index];
                    dmtx_.unlock();
                }                
                else{
                    data = tmp_;
                    return false;
                }
            }
    
            T& index(size_t i){
                
                if(i < size()){
                    dmtx_.lock();
                    int index = 0;
                    index = start_+i;
                    if(start_+i >= capacity_+1){
                        index = start_ + i - capacity_-1;
                    }
                    auto& tmp = datas_[index];
                    dmtx_.unlock();
                    return tmp;
                }
                else{
                    return tmp_;
                }
            }

            T& front(){
                if(empty()){
                    return tmp_;
                }
                return datas_[start_];
            }

            T& tail(){
                if(empty()){
                    return tmp_;
                }
                int index = end_ - 1;
                if(index == -1){
                    index = capacity_;
                }
                return datas_[index];
            }
    };
}

