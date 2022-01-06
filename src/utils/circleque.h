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

            bool push_front(){
                if(full()){
                    return false;
                }
                else{
                    --start_;
                    if(start_ == -1){
                        start_ = capacity_ - 1;
                    }
                }
            }

            void push_front_focus(){
                if(full()){
                    --end_;
                    --start_;
                    if(end_ == -1){
                        end_ = capacity_- 1;
                    }
                    if(start_ == -1){
                        start_ = capacity_ - 1;
                    }
                }
                else{
                    push_front();
                }
            }

            bool push_back(){
                if(full()){
                    return false;
                }
                else{
                    ++end_;
                    if(end_ == capacity_){
                        end_ = 0;
                    }
                }
            }

            void push_back_focus(){
                if(full()){
                    ++end_;
                    ++start_;
                    if(end_ == capacity_){
                        end_ = 0;
                    }
                    if(start_ == capacity_){
                        start_ = 0;
                    }
                }
                else{
                    push_back();
                }
            }

            bool pop_front(){
                if(empty()){
                    return false;
                }
                ++start;
                if(start_ == capacity_){
                    start_ = 0;
                }   
            }

            bool pop_back(){
                if(empty()){
                    return false;
                }
                --end_;
                if(end_ == -1){
                    end_ = capacity_-1;
                }
            }

            bool index(int i, T& data){
                if(i < 0 || i >= size()){
                    int index = 0;
                    index = start_+i;
                    if(start_+i >= capacity_){
                        index = start_ + i - capacity_;
                    }
                    data = datas_[index];
                }
                else{
                    data = T();
                    return false;
                }
            }
    };
}

