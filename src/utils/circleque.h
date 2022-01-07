#pragma once

#include <vector>
#include <atomic>
#include <mutex>
namespace svv_fusion{
    template <class T>
    class CircleQue{
        private:
            int capacity_;
            std::mutex dmtx_;
            std::vector<T> datas_;
            std::atomic_int start_;
            std::atomic_int end_;
        public:
            CircleQue(const int sz = 100):capacity_(sz){
                datas_.resize(capacity_);
                start_.store(0);
                end_.store(0);
            }
            int size(){
                if(empty())
                    return 0;
                if(full())
                    return capacity_;
                if(end_ >= start_)
                    return end_ - start_ + 1;
                else
                    return end_ + capacity_ - start_ + 1;
            }
            bool empty(){
                if(start_.load() == end_.load()){
                    return true;
                }
                return false;
            }

            bool full(){
                if(start_.load() == 0 || end_.load() == capacity_-1)
                    return true;
                if(start_.load() - 1 == end_.load() || end_.load() + 1 == start_.load())
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
                        start_ = capacity_ - 1;
                    }
                    datas_[start_] = data;
                }
                dmtx_.unlock();
            }

            void push_front_focus(T& data){
                dmtx_.lock();
                if(full()){
                    --end_;
                    --start_;
                    if(end_ == -1){
                        end_ = capacity_- 1;
                    }
                    if(start_ == -1){
                        start_ = capacity_ - 1;
                    }
                    datas_[start_] = data;
                }
                else{
                    push_front(data);
                }
                dmtx_.unlock();
            }

            bool push_back(T& data){
                if(full()){
                    return false;
                }
                else{
                    dmtx_.lock();
                    ++end_;
                    if(end_ == capacity_){
                        end_ = 0;
                    }
                    datas_[end_] = data;
                    dmtx_.unlock();
                }
                return true;
            }

            void push_back_focus(T& data){
                dmtx_.lock();
                if(full()){
                    ++end_;
                    ++start_;
                    if(end_ == capacity_){
                        end_ = 0;
                    }
                    if(start_ == capacity_){
                        start_ = 0;
                    }
                    datas_[end_] = data;
                }
                else{
                    push_back(data);
                }
                dmtx_.unlock();
            }

            bool pop_front(){
                if(empty()){
                    return false;
                }
                dmtx_.lock();

                ++start_;
                if (start_ == capacity_)
                {
                    start_ = 0;
                }
                dmtx_.unlock();
            }

            bool pop_back(){

                if(empty()){
                    return false;
                }
                dmtx_.lock();
                --end_;
                if(end_ == -1){
                    end_ = capacity_-1;
                }
                dmtx_.unlock();
            }

            bool index(size_t i, T& data){
                
                if(i < 0 || i >= size()){
                    dmtx_.lock();
                    int index = 0;
                    index = start_+i;
                    if(start_+i >= capacity_){
                        index = start_ + i - capacity_;
                    }
                    data = datas_[index];
                    dmtx_.unlock();
                }                
                else{
                    data = T();
                    return false;
                }
                
            }
    
            T index(size_t i){
                
                if(i < size()){
                    dmtx_.lock();
                    int index = 0;
                    index = start_+i;
                    if(start_+i >= capacity_){
                        index = start_ + i - capacity_;
                    }
                    auto& tmp = datas_[index];
                    dmtx_.unlock();
                    return tmp;
                    
                }
                else{
                    return T();
                }
            }
    };
}

