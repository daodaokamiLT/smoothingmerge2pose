#include "svvfusion.h"
#include "svvfactors.h"
#include <ceres/ceres.h>

namespace svv_fusion{
    VIOVPSFusion::VIOVPSFusion(const int& vioquesize, const int& vpsquesize, 
                    const int& matchsize):min_opt_size_(matchsize/2), 
                    vio_poses_(vioquesize), vps_poses_(vpsquesize), 
                    opt_vio_poses_(vioquesize), viovps_matches_(matchsize){
        initialized_.store(false);
        opt_thread_ = std::thread(&VIOVPSFusion::RunOptical, this);
    }

    VIOVPSFusion::~VIOVPSFusion(){
        opt_thread_.detach();
    }
    // 可能出现两次一样的match同时push 进入队列中！！！
    void VIOVPSFusion::PushVIOPose(Posed_t &pose) {
        bool moveflag = false;
        if(vio_poses_.full()){
            moveflag = true;
        }
        vio_poses_.push_back_focus(pose);        
        if(moveflag){
            MatchesFirstD1();
        }
        // todo if has new matches in this que...
        int start_index = 0;
        double time = pose.timestamp;
        if(!viovps_matches_.empty()){
            start_index = viovps_matches_.index(viovps_matches_.size()-1).second;   
            ++start_index;
        }
        int index = findTimeStampInVPS(start_index, time);
        if(index >= 0){
            newvps_match_.store(true);
            std::pair<size_t, size_t> match = std::make_pair(vio_poses_.size()-1, index);
            mlocker_.lock();
            if(match.first != viovps_matches_.tail().first){
                viovps_matches_.push_back_focus(match);
            }
            mlocker_.unlock();
        }
    }
    void VIOVPSFusion::PushVPSPose(Posed_t &pose) {
        bool moveflag = false;
        if(vps_poses_.full()){
            moveflag = true;
        }
        vps_poses_.push_back_focus(pose);
        if(moveflag){
            MatchesSecondD1();
        }
        // todo if has new matches in this que...
        int start_index = 0;
        double time = pose.timestamp;
        if(!viovps_matches_.empty()){
            start_index = viovps_matches_.index(viovps_matches_.size()-1).first;   
            ++start_index;
        }
        int index = findTimeStampInVIO(start_index, time);
        if(index >= 0){
            newvps_match_.store(true);
            std::pair<size_t, size_t> match = std::make_pair(index, vps_poses_.size()-1);
            mlocker_.lock();
            if(match.first != viovps_matches_.tail().first){
                viovps_matches_.push_back_focus(match);
            }
            mlocker_.unlock(); 
        }
    }

    bool VIOVPSFusion::SimpleInitializeByVPSPose()
    {
        // 假定，第一次不会偏的离谱 ？？？
        if(!viovps_matches_.empty()){
            initialized_.store(true);
            T_wvps_wvio_.block<3, 1>(0, 3) = Vec3d::Zero();
            T_wvps_wvio_.block<3, 3>(0, 0) = Mat3d::Identity();

            // @todo find the matches vps and vio
            auto viop = vio_poses_.index(viovps_matches_.index(0).first);
            auto vpsp = vps_poses_.index(viovps_matches_.index(0).second);

            T_wvps_wvio_.block<3,3>(0,0) = vpsp.q_wc.toRotationMatrix() * viop.q_wc.toRotationMatrix().transpose();
            T_wvps_wvio_.block<3,1>(0,3) = vpsp.q_wc.toRotationMatrix() * (-viop.q_wc.toRotationMatrix().transpose() * viop.t_wc) + vpsp.t_wc;

            initialized_.store(true);
        }
        return true;
    }

    bool VIOVPSFusion::SimpleInitializePose()
    {   
        // need at least 4 point and hasn't in one line
        // analysic solver
        Vec3d vio_center = Vec3d::Zero();
        Vec3d vps_center = Vec3d::Zero();
        std::vector<Posed_t> temp_vio_poses;
        std::vector<Posed_t> temp_vps_poses;
        {
            mlocker_.lock();
            int ms = viovps_matches_.size();
            if (ms < min_opt_size_)
            {
                return false;
            }
            temp_vio_poses.resize(ms);
            temp_vps_poses.resize(ms);
            // @todo 选择尽量远的几个pose 来做ICP
            for (int i = 0; i < ms; ++i)
            {
                temp_vio_poses[i] = vio_poses_.index(viovps_matches_.index(i).first);
                temp_vps_poses[i] = vps_poses_.index(viovps_matches_.index(i).second);
            }
            mlocker_.unlock();
        }

        int real_size =temp_vio_poses.size(); 
        // check in one line 
        Mat3d linecheck;
        // project to xoy
        bool isallinoneline = false;
        for(int i=0; i<real_size-2; ++i){
            Vec3d a = temp_vio_poses[i].t_wc;
            a[2] = 1;
            Vec3d b = temp_vio_poses[i+1].t_wc;
            b[2] = 1;
            Vec3d c = temp_vio_poses[i+2].t_wc;
            c[2] = 1;

            linecheck.block<1,3>(0,0) = a.transpose();
            linecheck.block<1,3>(1,0) = b.transpose();
            linecheck.block<1,3>(2,0) = c.transpose();

            if(std::fabs(linecheck.determinant()) > 1e-6){
                // is not in a line
                isallinoneline = true;
            }
        }
        // project to xoz
        if(!isallinoneline){
            for(int i=0; i<real_size-2; ++i){
                Vec3d a = temp_vio_poses[i].t_wc;
                a[1] = a[2]; a[2] = 1;
                Vec3d b = temp_vio_poses[i + 1].t_wc;
                b[1] = b[2]; b[2] = 1;
                Vec3d c = temp_vio_poses[i + 2].t_wc;
                c[2] = c[2]; c[2] = 1;

                linecheck.block<1, 3>(0, 0) = a.transpose();
                linecheck.block<1, 3>(1, 0) = b.transpose();
                linecheck.block<1, 3>(2, 0) = c.transpose();

                if (std::fabs(linecheck.determinant()) > 1e-6)
                {
                    // is not in a line
                    isallinoneline = true;
                }
            }
        }
        // project to yoz
        if(!isallinoneline){
            for(int i=0; i<real_size-2; ++i){
                Vec3d a = temp_vio_poses[i].t_wc;
                a[0] = a[1]; a[1] = a[2]; a[2] = 1;
                Vec3d b = temp_vio_poses[i + 1].t_wc;
                b[0] = b[1]; b[1] = b[2]; b[2] = 1;
                Vec3d c = temp_vio_poses[i + 2].t_wc;
                c[0] = c[1]; c[2] = c[2]; c[2] = 1;

                linecheck.block<1, 3>(0, 0) = a.transpose();
                linecheck.block<1, 3>(1, 0) = b.transpose();
                linecheck.block<1, 3>(2, 0) = c.transpose();

                if (std::fabs(linecheck.determinant()) > 1e-6)
                {
                    // is not in a line
                    isallinoneline = true;
                }
            }
        }

        if(!isallinoneline){
            return false;
        }
        for(int i=0; i<real_size; ++i){
            vio_center += temp_vio_poses[i].t_wc;
            vps_center += temp_vps_poses[i].t_wc;
        }

        vio_center /= real_size;
        vps_center /= real_size;

        for(int i=0; i<real_size; ++i){
            temp_vio_poses[i].t_wc = temp_vio_poses[i].t_wc - vio_center;
            temp_vps_poses[i].t_wc = temp_vps_poses[i].t_wc - vps_center;
        }

        // vps to vio 
        // the function: i is vps, j is vio
        Mat3d H = Mat3d::Zero();
        for(int i=0; i<real_size; ++i){
            H += temp_vps_poses[i].t_wc * temp_vio_poses[i].t_wc.transpose();
        }
        Eigen::JacobiSVD<Mat3d> svd(H, Eigen::ComputeFullU|Eigen::ComputeFullV);
        Mat3d U = svd.matrixU();
        Mat3d V = svd.matrixV();

        T_wvps_wvio_.block<3,3>(0,0) = (U * V.transpose()).transpose();
        T_wvps_wvio_.block<3,1>(0,3) = vio_center - T_wvps_wvio_.block<3,3>(0,0).transpose() * vps_center;
        std::cout<<"T_vps_vio\n "<<T_wvps_wvio_<<std::endl;
        
        initialized_.store(true);
        return true;
    }

    bool VIOVPSFusion::RansacInitializePose() {
        printf("the ransac method havn't completed...\n");
        return false;
    }

    void VIOVPSFusion::RunOptical() {
        CircleQue<Posed_t> tmp_vioposes(viovps_matches_.capacity());
        CircleQue<Posed_t> tmp_vpsposes(viovps_matches_.capacity());
        while(true){
            printf("run in Optical....\n");
            if(newvps_match_.load()){
                newvps_match_.store(false);
                mlocker_.lock();
                std::pair<size_t, size_t> pair = viovps_matches_.tail();
                tmp_vioposes.push_back_focus(vio_poses_.index(pair.first));
                tmp_vpsposes.push_back_focus(vps_poses_.index(pair.second));
                mlocker_.unlock();
                Posed_t last_wvps_viopose, wvps_viopose;
                last_wvps_viopose.t_wc;
                last_wvps_viopose.q_wc;

                wvps_viopose.t_wc;
                wvps_viopose.q_wc;
                // T_wvps_wvio  should also be opticated
                Optical(tmp_vioposes, tmp_vpsposes, last_wvps_viopose, wvps_viopose);
            }
            else{
                std::chrono::milliseconds dura(2000);
                std::this_thread::sleep_for(dura);
            }
        }
    }

    void VIOVPSFusion::Optical(CircleQue<Posed_t>& vioposes, CircleQue<Posed_t>& vpsposes, 
                               Posed_t& last_pose, Posed_t& cur_pose){
        // a another element should also be opti ed !!! T_wvps_wvio
        // the last_pose always the tail()-1, cur_pose is the tail pose
        // 优化的变量实际是 delta_pose
        ceres::Problem problem;
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
        options.max_num_iterations = 5;
        ceres::Solver::Summary summary;
        ceres::LossFunction *loss_function;
        loss_function = new ceres::HuberLoss(1.0);
        ceres::LocalParameterization* local_parameteriztion = new ceres::EigenQuaternionParameterization();
        // fix front 4/5 opt the tail 1/5 pose
        // always change the T_wvps_wvio
        
        int deppart = min_opt_size_*4.0/5.0;
        if(vioposes.size() < deppart){
            // when pose not enough ... 
            return;
        }
        
        Quaterniond qvec_wvps_wvio (T_wvps_wvio_.block<3,3>(0,0));
        Vec3d t_wvps_wvio = T_wvps_wvio_.block<3,1>(0,3);

        problem.AddParameterBlock(qvec_wvps_wvio.coeffs().data(), 4, local_parameteriztion);
        problem.AddParameterBlock(t_wvps_wvio.data(), 3);

        for(int i=0; i<deppart; ++i){
            auto& viop = vioposes.index(i);
            auto& vpsp = vioposes.index(i);
            if(!viop.t_wc.hasNaN() && !vpsp.t_wc.hasNaN()){
                problem.AddParameterBlock(viop.q_wc.coeffs().data(), 4, local_parameteriztion);
                problem.AddParameterBlock(viop.t_wc.data(), 3);

                problem.SetParameterBlockConstant(viop.q_wc.coeffs().data());
                problem.SetParameterBlockConstant(viop.t_wc.data());
                
                // vps result 是可以优化一下的！！！
                problem.AddParameterBlock(vpsp.q_wc.coeffs().data(), 4, local_parameteriztion);
                problem.AddParameterBlock(vpsp.t_wc.data(), 3);
            }
            else{
                printf("error before depart data error...\n");
                return;
            }
        }
        for(int i=deppart; i<vioposes.size(); ++i){
            auto& viop = vioposes.index(i);
            auto& vpsp = vpsposes.index(i);

            if(!viop.t_wc.hasNaN() && !vpsp.t_wc.hasNaN()){
                problem.AddParameterBlock(viop.q_wc.coeffs().data(), 4, local_parameteriztion);
                problem.AddParameterBlock(viop.t_wc.data(), 3);

                problem.AddParameterBlock(vpsp.q_wc.coeffs().data(), 4, local_parameteriztion);
                problem.AddParameterBlock(vpsp.t_wc.data(), 3);
            }
            else{
                printf("error after depart data error...\n");
                return;
            }
        }

        auto &viop0 = vioposes.index(0);
        auto &vpsp0 = vpsposes.index(0);

        ceres::CostFunction* cviovps_function0 = RelativeCRTError::Create(vpsp0.t_wc[0], vpsp0.t_wc[1], vpsp0.t_wc[2], 
                            vpsp0.q_wc.w(), vpsp0.q_wc.x(), vpsp0.q_wc.y(), vpsp0.q_wc.z(), 
                            0.1, 0.01);

        problem.AddResidualBlock(cviovps_function0, loss_function, t_wvps_wvio.data(), qvec_wvps_wvio.coeffs().data(), 
                                viop0.t_wc.data(), viop0.q_wc.coeffs().data());

        for(int i=1; i<vioposes.size(); ++i){
            auto& viop1 = vioposes.index(i);
            auto& vpsp1 = vpsposes.index(i);

            ceres::CostFunction* cviovps_function1 = RelativeCRTError::Create(vpsp1.t_wc[0], vpsp1.t_wc[1], vpsp1.t_wc[2], 
                            vpsp1.q_wc.w(), vpsp1.q_wc.x(), vpsp1.q_wc.y(), vpsp1.q_wc.z(), 
                            0.1, 0.01);

            problem.AddResidualBlock(cviovps_function1, loss_function, t_wvps_wvio.data(), qvec_wvps_wvio.coeffs().data(), 
                                viop1.t_wc.data(), viop1.q_wc.coeffs().data());

            Posed_t deltavps01;
            deltaPosed(vpsp0, vpsp1, deltavps01);
            ceres::CostFunction* viovps_function = RelativeRTError::Create(deltavps01.t_wc[0], deltavps01.t_wc[1], deltavps01.t_wc[2],
                                                        deltavps01.q_wc.w(), deltavps01.q_wc.x(), deltavps01.q_wc.y(), deltavps01.q_wc.z(),
                                                        0.1, 0.01);
            problem.AddResidualBlock(viovps_function, NULL, viop0.q_wc.coeffs().data(), viop0.t_wc.data(), 
                                viop1.q_wc.coeffs().data(), viop1.t_wc.data());

            viop0 = viop1;
            vpsp0 = vpsp1;
        }
        ceres::Solve(options, &problem, &summary);
        
        T_wvps_wvio_.block<3,3>(0,0) = qvec_wvps_wvio.toRotationMatrix();
        T_wvps_wvio_.block<3,1>(0,3) = t_wvps_wvio;

        opt_vio_poses_.push_back_focus(vioposes.tail());
    }

    void VIOVPSFusion::GetWVPS_VIOPose(Posed_t& pose) {
        pose.timestamp = cur_timestamp_;
        pose.t_wc = t_wvps_vio_;
        pose.q_wc = q_wvps_vio_;
    }
    

    int VIOVPSFusion::findTimeStampInVPS(int start_index, double timestamp){
        if(start_index>=0 && start_index < vps_poses_.size()){
            for(size_t i=start_index; i<vps_poses_.size(); ++i){
                if(vps_poses_.index(i).timestamp == timestamp){
                    return i;
                }
            }
        }
        return -1;
    }
    int VIOVPSFusion::findTimeStampInVIO(int start_index, double timestamp){
        if(start_index>=0 && start_index < vio_poses_.size()){
            for(size_t i=start_index; i<vio_poses_.size(); ++i){
                if(vio_poses_.index(i).timestamp == timestamp)
                    return i;
            }
        }
        return -1;
    }

    void VIOVPSFusion::MatchesFirstD1(){
        mlocker_.lock();
        int start = 0;
        if (viovps_matches_.index(start).first == 0)
        {
            viovps_matches_.pop_front();
            start = 1;
        }
        for (int i = start; i < viovps_matches_.size(); ++i)
        {
            viovps_matches_.index(i).first -= 1;
        }
        mlocker_.unlock();
    }
    void VIOVPSFusion::MatchesSecondD1(){
        mlocker_.lock();
        int start = 0;
        if(viovps_matches_.index(start).second == 0){
            viovps_matches_.pop_front();
            start = 1;
        }
        for(int i=start; i<viovps_matches_.size(); ++i){
            viovps_matches_.index(i).second -= 1;
        }
        mlocker_.unlock();
    }
}