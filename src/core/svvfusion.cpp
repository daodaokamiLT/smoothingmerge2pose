#include "svvfusion.h"

namespace svv_fusion{
    VIOVPSFusion::VIOVPSFusion(const int& vioquesize, const int& vpsquesize, 
                    const int& matchsize):min_opt_size_(matchesize), 
                    vio_poses_(vioquesize), vps_poses_(vpsquesize), 
                    matches_(matchsize){
        
        
        
    }

    // 初始化是否需要20个连续的keyframe进行match
    void VIOVPSFusion::SetVIOVPS_Initialize(int size) {


    }

    void VIOVPSFusion::PushVIOPose(Posed_t &pose) {



    }
    void VIOVPSFusion::PushVPSPose(Posed_t &pose) {



    }

    bool VIOVPSFusion::InitializePose() {



    }

    bool VIOVPSFusion::RansacInitializePose() {


    }

    void VIOVPSFusion::RunOptical() {


    }
    
    void VIOVPSFusion::GetWVPS_VIOPose() {



    }
}