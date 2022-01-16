#include "svvfusion.h"
#include "svvfactors.h"
#include <ceres/ceres.h>
#include <fstream>

namespace svv_fusion{
    VIOVPSFusion::VIOVPSFusion(const int& vioquesize, const int& vpsquesize, 
                    const int& matchsize):min_opt_size_(matchsize/2), 
                    vio_poses_(vioquesize), vps_poses_(vpsquesize), 
                    opt_vio_poses_(vioquesize), viovps_matches_(matchsize){
        matches_index_ = 0;
        newvps_match_.store(false);
        initialized_.store(false);
        opt_thread_ = std::thread(&VIOVPSFusion::RunOptical, this);
    }

    VIOVPSFusion::~VIOVPSFusion(){
        opt_thread_.detach();
        
        std::ofstream foutC("/home/lut/Desktop/evo/optvio.csv", std::ios::app);
        // foutC.setf(std::ios::fixed, std::ios::floatfield);
        //     foutC.precision(0);
        //     foutC << optpose.timestamp << ",";
        //     foutC.precision(5);
        //     foutC << optpose.t_wc[0] << ","
        //           << optpose.t_wc[1] << ","
        //           << optpose.t_wc[2] << ","
        //           << optpose.q_wc.w() << ","
        //           << optpose.q_wc.x() << ","
        //           << optpose.q_wc.y() << ","
        //           << optpose.q_wc.z() << std::endl;
        foutC.close();
    }
    // 可能出现两次一样的match同时push 进入队列中！！！
    // 可能出现两次一样的match同时push 进入队列中！！！ // 这种find matches 方法太差了！！！！
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
        int viostart = 0;
        int vpsstart = 0;
        if(!viovps_matches_.empty()){
            auto m = viovps_matches_.index(matches_index_);
            viostart = m.first+1;
            vpsstart = m.second+1;    
        }
        while(viostart < vio_poses_.size() && vpsstart < vps_poses_.size()){
            Posed_t& viop = vio_poses_.index(viostart);
            Posed_t& vpsp = vps_poses_.index(vpsstart);
            if(viop.timestamp == vpsp.timestamp){
                if(viostart == -1){
                    printf("error, startvio index is -1...\n");
                    exit(-1);
                }
                std::pair<size_t, size_t> match = std::make_pair(viostart, vpsstart);
                viovps_matches_.push_back_focus(match);
                matches_index_ = viovps_matches_.size()-1;
                newvps_match_.store(true);
                break;
            }
            else if(viop.timestamp < vpsp.timestamp){
                ++viostart;
            }
            else if(viop.timestamp > vpsp.timestamp){
                ++vpsstart;
            }
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
        int viostart = 0;
        int vpsstart = 0;
        if(!viovps_matches_.empty()){
            auto m = viovps_matches_.index(matches_index_);
            viostart = m.first+1;
            vpsstart = m.second+1;       
        }
        while(viostart < vio_poses_.size() && vpsstart < vps_poses_.size()){
            Posed_t& viop = vio_poses_.index(viostart);
            Posed_t& vpsp = vps_poses_.index(vpsstart);
            if(viop.timestamp == vpsp.timestamp){
                std::pair<size_t, size_t> match = std::make_pair(viostart, vpsstart);
                viovps_matches_.push_back_focus(match);
                matches_index_ = viovps_matches_.size()-1;
                newvps_match_.store(true);
                break;
            }
            else if(viop.timestamp < vpsp.timestamp){
                ++viostart;
            }
            else if(viop.timestamp > vpsp.timestamp){
                ++vpsstart;
            }
        }
    }

    bool VIOVPSFusion::SimpleInitializeByVPSPose()
    {
        // 假定，第一次不会偏的离谱 ？？？
        if(!viovps_matches_.empty()){
            if(viovps_matches_.size() == 1){
                initialized_.store(true);
                T_wvps_wvio_.block<3, 1>(0, 3) = Vec3d::Zero();
                T_wvps_wvio_.block<3, 3>(0, 0) = Mat3d::Identity();

                // @todo find the matches vps and vio
                auto viop = vio_poses_.index(viovps_matches_.index(0).first);
                auto vpsp = vps_poses_.index(viovps_matches_.index(0).second);

                T_wvps_wvio_.block<3, 3>(0, 0) = vpsp.q_wc.toRotationMatrix() * viop.q_wc.toRotationMatrix().transpose();
                T_wvps_wvio_.block<3, 1>(0, 3) = vpsp.q_wc.toRotationMatrix() * (-viop.q_wc.toRotationMatrix().transpose() * viop.t_wc) + vpsp.t_wc;
            }
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
            if(newvps_match_.load()){
                printf("run in Optical matches size is %d....\n", viovps_matches_.size());
                newvps_match_.store(false);
                mlocker_.lock();
                std::pair<size_t, size_t> pair = viovps_matches_.tail();
                printf("pair.first second %d %d...\n", pair.first, pair.second);
                Posed_t viop = vio_poses_.index(pair.first);
                Posed_t vpsp = vps_poses_.index(pair.second);
                printf("match size %d, the match element ? ... %d %d %lf.\n",viovps_matches_.size(), 
                                    pair.first, pair.second, viop.timestamp);
                double t = viop.timestamp;
                if(t != vpsp.timestamp){
                    printf("error, match size %d, the match element didn't matched... %d %d|| %lf %lf\n", 
                                    viovps_matches_.size(), 
                                    pair.first, pair.second, 
                                    t, vpsp.timestamp);
                    exit(-1);
                }
                printf("tmp_vio and vps poses size is %d %d...\n", tmp_vioposes.size(), tmp_vpsposes.size());
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
                std::chrono::milliseconds dura(10);
                std::this_thread::sleep_for(dura);
            }
        }
    }

    bool SimpleAligned(CircleQue<Posed_t> &vp0, CircleQue<Posed_t> &vp1, Mat4d &T_01)
    {
        if(vp0.size() != vp1.size()){
            printf("aligned error....\n");
            return false;
        }
        // move 16cm, cannot make initialized?
        // 这应该是能够计算出deltapose的？？？
        auto deltavp0 = vp0.tail().t_wc - vp0.front().t_wc;
        std::cout<<"vps deltapose "<< deltavp0.transpose() <<" "<<deltavp0.norm()<<std::endl; 

        auto deltavp1 = vp1.tail().t_wc - vp1.front().t_wc;
        std::cout<<"vio deltapose "<< deltavp1.transpose() <<" "<<deltavp1.norm()<<std::endl; 
        // all element is matched
        T_01 = Mat4d::Identity();
        ceres::Problem problem;
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
        options.max_num_iterations = 5;
        ceres::Solver::Summary summary;
        ceres::LossFunction *loss_function;
        loss_function = new ceres::HuberLoss(1.0);
        ceres::LocalParameterization *local_parameteriztion = new ceres::EigenQuaternionParameterization();
        Quaterniond qvec_01 = Quaterniond::Identity();
        Vec3d t_01 = Vec3d::Zero();
        qvec_01 = vp0.index(0).q_wc.inverse() * vp1.index(0).q_wc;
        t_01 = vp0.index(0).q_wc.inverse() * (vp1.index(0).t_wc - vp0.index(0).t_wc);
        problem.AddParameterBlock(qvec_01.coeffs().data(), 4, local_parameteriztion);
        problem.AddParameterBlock(t_01.data(), 3);
        std::cout<<"vio vps pose size is ..."<<vp0.size()<<std::endl;
        for(int i=0; i<vp0.size(); ++i){
            auto p0 = vp0.index(i);
            auto p1 = vp1.index(i);

            std::cout<<i<<"viop: "<<p0<<std::endl;
            std::cout<<i<<"vpsp: "<<p1<<std::endl;
            ceres::CostFunction* init_function = ChangeCoordinateError::Create(p0.t_wc[0], p0.t_wc[1], p0.t_wc[2], 
            p0.q_wc.w(), p0.q_wc.x(), p0.q_wc.y(), p0.q_wc.z(), 1, 0.1,
            p1.t_wc[0], p1.t_wc[1], p1.t_wc[2], 
            p1.q_wc.w(), p1.q_wc.x(), p1.q_wc.y(), p1.q_wc.z(), 1, 0.1);

            problem.AddResidualBlock(init_function, NULL, qvec_01.coeffs().data(), t_01.data());
        }

        ceres::Solve(options, &problem, &summary);
        std::cout<<summary.BriefReport()<<std::endl;
        T_01.block<3,3>(0,0) = qvec_01.toRotationMatrix();
        T_01.block<3,1>(0,3) = t_01;

        auto p0 = vp0.tail();
        auto p1 = vp1.tail();
        auto delta_R = p0.q_wc.toRotationMatrix().transpose() * p1.q_wc.toRotationMatrix();
        auto delta_t = p0.q_wc.toRotationMatrix().transpose() * (p1.t_wc - p0.t_wc);

        std::cout<<T_01<<std::endl;
        std::cout<<delta_R <<"\n"<<delta_t<<std::endl;

        p0 = vp0.front();
        p1 = vp1.front();
        auto delta_R1 = p0.q_wc.toRotationMatrix().transpose() * p1.q_wc.toRotationMatrix();
        auto delta_t1 = p0.q_wc.toRotationMatrix().transpose() * (p1.t_wc - p0.t_wc);
        std::cout<<delta_R1 <<"\n"<<delta_t1<<std::endl;

        exit(-1);
        return true;
    }

    void VIOVPSFusion::Optical(CircleQue<Posed_t>& vioposes, CircleQue<Posed_t>& vpsposes, 
                               Posed_t& last_pose, Posed_t& cur_pose){
        // a another element should also be opti ed !!! T_wvps_wvio
        // the last_pose always the tail()-1, cur_pose is the tail pose
        // 优化的变量实际是 delta_pose
        printf("real run in Optical solve the problem.\n");
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
        int deppart = std::max(min_opt_size_, vioposes.size())*4.0/5.0;
        if(vioposes.size() < deppart){
            // when pose not enough ... 
            return;
        }
        else if(vioposes.size() == deppart){
            // cal the init coodinate change
            SimpleAligned(vpsposes, vioposes, T_wvps_wvio_);
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
                std::cout<<i<<"viop:"<<viop<<std::endl;
                std::cout<<i<<"vpsp:"<<vpsp<<std::endl;
                exit(0);
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
                std::cout<<i<<"viop:"<<viop<<std::endl;
                std::cout<<i<<"vpsp:"<<vpsp<<std::endl;
                exit(-1);
                return;
            }
        }

        // RelativeRT error + Rt error ?
        auto &viop0 = vioposes.index(0);
        auto &vpsp0 = vpsposes.index(0);
        ceres::CostFunction* vps_function = RTError::Create(vpsp0.t_wc[0], vpsp0.t_wc[1], vpsp0.t_wc[2], 
                                                        vpsp0.q_wc.w(), vpsp0.q_wc.x(), vpsp0.q_wc.y(), vpsp0.q_wc.z(), 0.1, 0.01);
        problem.AddResidualBlock(vps_function, loss_function, 
                            viop0.q_wc.coeffs().data(), viop0.t_wc.data());
        for(int i=1; i<vioposes.size(); ++i){
            auto& viop1 = vioposes.index(i);
            auto& vpsp1 = vpsposes.index(i);

            Mat4d wTi = Mat4d::Identity();
            Mat4d wTj = Mat4d::Identity();

            wTi.block<3, 3>(0, 0) = viop0.q_wc.toRotationMatrix();
            wTi.block<3, 1>(0, 3) = viop0.t_wc;

            wTj.block<3, 3>(0, 0) = viop1.q_wc.toRotationMatrix();
            wTj.block<3, 1>(0, 3) = viop1.t_wc;

            Mat4d iTj = wTi.inverse() * wTj;
            Quaterniond iQj(iTj.block<3, 3>(0, 0));
            Vec3d iPj = iTj.block<3, 1>(0, 3);

            ceres::CostFunction *vio_function = RelativeRTError::Create(iPj.x(), iPj.y(), iPj.z(),
                                                                        iQj.w(), iQj.x(), iQj.y(), iQj.z(), 0.1, 0.01);
            problem.AddResidualBlock(vio_function, NULL, vpsp0.q_wc.coeffs().data(), vpsp0.t_wc.data(),
                                     vpsp1.q_wc.coeffs().data(), vpsp1.t_wc.data());


            ceres::CostFunction* vps_function = RTError::Create(vpsp1.t_wc[0], vpsp1.t_wc[1], vpsp1.t_wc[2], 
                                                        vpsp1.q_wc.w(), vpsp1.q_wc.x(), vpsp1.q_wc.y(), vpsp1.q_wc.z(), 0.1, 0.01);
            problem.AddResidualBlock(vps_function, loss_function, 
                            viop1.q_wc.coeffs().data(), viop1.t_wc.data());                                                        
            viop0 = viop1;
            vpsp0 = vpsp1;
        }
        ceres::Solve(options, &problem, &summary);
        T_wvps_wvio_.block<3,3>(0,0) = qvec_wvps_wvio.toRotationMatrix();
        T_wvps_wvio_.block<3,1>(0,3) = t_wvps_wvio; 
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