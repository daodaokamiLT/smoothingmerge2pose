#pragma once

#include "../utils/basic_types.h"
#include "../utils/circleque.h"
#include <thread>
#include <mutex>
#include <atomic>
namespace svv_fusion{
    class VIOVPSFusion{
        private:
            int min_opt_size_;
            AlignedMap<double, Posed_t> allopted_vio_poses_; 

            CircleQue<Posed_t> vio_poses_;
            CircleQue<Posed_t> vps_poses_;
            CircleQue<Posed_t> opt_vio_poses_;

            int findTimeStampInVPS(int start_index, double timestamp);
            int findTimeStampInVIO(int start_index, double timestamp);
            std::atomic_bool newvps_match_;

            void MatchesFirstD1();
            void MatchesSecondD1();
            double last_timestamp_, cur_timestamp_;
            Vec3d t_wvps_vio_;
            Vec3d t_wvps_vio_last_;
            Quaterniond q_wvps_vio_;
            Quaterniond q_wvps_vio_last_; 
            // for coordinate change !!!
            Mat4d T_wvps_wvio_ = Mat4d::Identity();
            std::mutex mlocker_;
            std::atomic_int matches_index_;
            CircleQue<std::pair<size_t, size_t>> viovps_matches_;
            std::thread opt_thread_;
            std::atomic_bool initialized_;
        public:
            double last_vio_timestamp_;
            double last_vps_timestamp_;
            double delta_timestamp_vio2vps_;

            VIOVPSFusion(const int& vioquesize = 100, const int& vpsquesize = 50, 
                    const int& matchsize = 50);
            ~VIOVPSFusion();
            // 初始化是否需要20个连续的keyframe进行match
            void PushVIOPose(Posed_t& pose);
            void PushVPSPose(Posed_t& pose);

            bool SimpleInitializeByVPSPose();
            bool SimpleInitializePose();
            bool RansacInitializePose();

            void RunOptical();
            void Optical(CircleQue<Posed_t>& vioposes, CircleQue<Posed_t>& vpsposes, 
                               Posed_t& last_pose, Posed_t& cur_pose);

            bool IsInitialization(){
                return initialized_.load();
            }
    };
}