#pragma once

#include "../utils/basic_types.h"
#include "../utils/circleque.h"

namespace svv_fusion{
    class VIOVPSFusion{
        private:
            int min_opt_size_;
            CircleQue<Posed_t> vio_poses_;
            CircleQue<Posed_t> vps_poses_;
            Vec3d t_wvps_vio_;
            Quaterniond wvps_vio_;
            // for coordinate change !!!
            Mat4d T_wvps_wvio_;
            CircleQue<std::pair<size_t, size_t>> matches_;

        public:
            double last_vio_timestamp_;
            double last_vps_timestamp_;
            double delta_timestamp_vio2vps_;

            VIOVPSFusion(const int& vioquesize = 100, const int& vpsquesize = 50, 
                    const int& matchsize = 50);

            // 初始化是否需要20个连续的keyframe进行match
            void SetVIOVPS_Initialize(int size = 20);
            void PushVIOPose(Posed_t& pose);
            void PushVPSPose(Posed_t& pose);

            bool InitializePose();
            bool RansacInitializePose();

            void RunOptical();
            void GetWVPS_VIOPose();
    };
}