#include "vinsfusion.h"
#include "vinsfactors.h"
#include <fstream>
#include <unistd.h>
namespace svv_fusion{
namespace vins_fusion{
    VIOVPSFusion2::VIOVPSFusion2(const int& vioquesize, const int& vpsquesize, 
                    const int& matchsize):min_opt_size_(matchsize/2), 
                    vio_poses_(vioquesize), vps_poses_(vpsquesize), 
                    opt_vio_poses_(vioquesize), viovps_matches_(matchsize){\
        matches_index_ = 0;
        initialized_.store(false);
        newvps_match_.store(false);
        opt_thread_ = std::thread(&VIOVPSFusion2::RunOptical, this);
    }

    VIOVPSFusion2::~VIOVPSFusion2(){
        opt_thread_.detach();
        std::ofstream foutC("/home/lut/Desktop/evo/optvinsvio.csv", std::ios::app);
        foutC.setf(std::ios::fixed, std::ios::floatfield);
        if(foutC.is_open()){
            for (auto iter : wvps_vioposeMap_)
            {
                if (iter.second.valued())
                {
                    foutC.precision(0);
                    foutC << iter.first << ",";
                    foutC.precision(5);
                    foutC << iter.second.t_wc[0] << ","
                          << iter.second.t_wc[1] << ","
                          << iter.second.t_wc[2] << ","
                          << iter.second.q_wc.w() << ","
                          << iter.second.q_wc.x() << ","
                          << iter.second.q_wc.y() << ","
                          << iter.second.q_wc.z() << std::endl;
                }
            }
            foutC.close();
        }
    }
    // 可能出现两次一样的match同时push 进入队列中！！！ // 这种find matches 方法太差了！！！！
    void VIOVPSFusion2::PushVIOPose(Posed_t &pose) {
        bool moveflag = false;
        if(vio_poses_.full()){
            moveflag = true;
        }
        vio_poses_.push_back_focus(pose);        
        if(moveflag){
            printf("move flag....\n");
            MatchesFirstD1();
        }
        // todo if has new matches in this que...
        int viostart = 0;
        int vpsstart = 0;
        
        if(!viovps_matches_.empty()){
            mlocker_.lock();
            auto m = viovps_matches_.index(matches_index_);
            mlocker_.unlock();
            viostart = m.first+1;
            vpsstart = m.second+1;    
        }
        while(viostart < vio_poses_.size() && vpsstart < vps_poses_.size()){
            Posed_t& viop = vio_poses_.index(viostart);
            Posed_t& vpsp = vps_poses_.index(vpsstart);
            if(viop.timestamp == vpsp.timestamp){
                mlocker_.lock();
                std::pair<size_t, size_t> match = std::make_pair(viostart, vpsstart);
                viovps_matches_.push_back_focus(match);
                matches_index_ = viovps_matches_.size()-1;
                mlocker_.unlock();
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
    void VIOVPSFusion2::PushVPSPose(Posed_t &pose) {
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
                mlocker_.lock();
                std::pair<size_t, size_t> match = std::make_pair(viostart, vpsstart);
                viovps_matches_.push_back_focus(match);
                mlocker_.unlock();

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

    bool VIOVPSFusion2::SimpleInitializeByVPSPose()
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

    bool VIOVPSFusion2::SimpleInitializePose()
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

    bool VIOVPSFusion2::RansacInitializePose() {
        printf("the ransac method havn't completed...\n");
        return false;
    }

    void VIOVPSFusion2::RunOptical() {
        while(true){
            if(newvps_match_.load()){
                newvps_match_.store(false);
                auto m = viovps_matches_.tail();
                Posed_t viop = vio_poses_.index(m.first);
                Posed_t vpsp = vps_poses_.index(m.second);
                printf("match size %d, the match element ? ... %d %d %lf.\n",viovps_matches_.size(), 
                                    m.first, m.second, viop.timestamp);
                double t = viop.timestamp;
                if(t != vpsp.timestamp){
                    
                    printf("error, match size %d, the match element didn't matched... %d %d|| %lf %lf\n", 
                                    viovps_matches_.size(), 
                                    m.first, m.second, 
                                    t, vpsp.timestamp);
                    exit(-1);
                }
                vioposeMap_[t] = viop;
                vpsposeMap_[t] = vpsp;
                Posed_t wvpsvioPose;
                ChangeCoordinate(T_wvps_wvio_, viop, wvpsvioPose);
                wvps_vioposeMap_[t] = wvpsvioPose;
                Optical();
            }
            else{
                std::chrono::milliseconds dura(10);
                std::this_thread::sleep_for(dura);
            }
        }
    }

    void VIOVPSFusion2::Optical(){
        ceres::Problem problem;
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
        options.max_num_iterations = 5;
        ceres::Solver::Summary summary;
        ceres::LossFunction* loss_function;
        loss_function = new ceres::HuberLoss(1.0);
        ceres::LocalParameterization* local_parameterization = new ceres::EigenQuaternionParameterization();

        int length = vioposeMap_.size();
        AlignedMap<double, Posed_t>::iterator iter;
        iter = wvps_vioposeMap_.begin();
        for(int i=0; i<length; i++, iter++){
            problem.AddParameterBlock(iter->second.q_wc.coeffs().data(), 4, local_parameterization);
            problem.AddParameterBlock(iter->second.t_wc.data(), 3);
        }

        AlignedMap<double, Posed_t>::iterator iterVIO, iterVIONext, iterVPS, iterwvpsVIO, iterwvpsVIONext;
        int i=0; 
        double last_timestamp = 0;
        for(iterVIO = vioposeMap_.begin(), iterwvpsVIO = wvps_vioposeMap_.begin(); 
                        iterVIO!=vioposeMap_.end(); 
                        iterVIO++, iterwvpsVIO++, i++){
            iterVIONext = iterVIO;
            iterVIONext++;

            iterwvpsVIONext = iterwvpsVIO;
            iterwvpsVIONext++;
            if(iterVIONext != vioposeMap_.end()){
                Mat4d wTi = Mat4d::Identity();
                Mat4d wTj = Mat4d::Identity();

                wTi.block<3,3>(0,0) = iterVIO->second.q_wc.toRotationMatrix();
                wTi.block<3,1>(0,3) = iterVIO->second.t_wc;

                wTj.block<3,3>(0,0) = iterVIONext->second.q_wc.toRotationMatrix();
                wTj.block<3,1>(0,3) = iterVIONext->second.t_wc;

                Mat4d iTj = wTi.inverse() * wTj;
                Quaterniond iQj(iTj.block<3,3>(0,0));
                Vec3d iPj = iTj.block<3,1>(0,3);

                ceres::CostFunction* vio_function = RelativeRTError::Create(iPj.x(), iPj.y(), iPj.z(), 
                                                                    iQj.w(), iQj.x(), iQj.y(), iQj.z(), 0.1, 0.01);
                problem.AddResidualBlock(vio_function, NULL, iterwvpsVIO->second.q_wc.coeffs().data(), iterwvpsVIO->second.t_wc.data(), 
                                        iterwvpsVIONext->second.q_wc.coeffs().data(), iterwvpsVIONext->second.t_wc.data());
            }

            double t = iterVIO->first;
            iterVPS = vpsposeMap_.find(t);
            if(iterVPS != vpsposeMap_.end()){
                ceres::CostFunction* vps_function = RTError::Create(iterVPS->second.t_wc[0], iterVPS->second.t_wc[1], iterVPS->second.t_wc[2], 
                                                        iterVPS->second.q_wc.w(), iterVPS->second.q_wc.x(), iterVPS->second.q_wc.y(), iterVPS->second.q_wc.z(), 0.1, 0.01);
                problem.AddResidualBlock(vps_function, loss_function, iterwvpsVIO->second.q_wc.coeffs().data(), iterwvpsVIO->second.t_wc.data());
            }
        }

        ceres::Solve(options, &problem, &summary);
        // 直接就是整个map的优化了
        // 这里实时的更新一下wvps_wvio 的变化
        iter = wvps_vioposeMap_.end();
        iter--;
        auto itervio = vioposeMap_.end();
        itervio--;
        if(iter->first != itervio->first){
            printf("error, timestamp not eqauls...\n");
            exit(-1);
        } 
        Mat4d T_wvps_vio = Mat4d::Identity();
        Mat4d T_wvio_vio = Mat4d::Identity();
        T_wvps_vio.block<3,3>(0,0) = iter->second.q_wc.toRotationMatrix();
        T_wvps_vio.block<3,1>(0,3) = iter->second.t_wc;
        T_wvio_vio.block<3,3>(0,0) = itervio->second.q_wc.toRotationMatrix();
        T_wvio_vio.block<3,1>(0,3) = itervio->second.t_wc;
        T_wvps_wvio_ = T_wvps_vio * T_wvio_vio.inverse();
        // std::cout<<"the T_wvps_wvio change method is\n"<<T_wvps_wvio_<<std::endl;
    }
    
    int VIOVPSFusion2::findTimeStampInVPS(int start_index, double timestamp){

        if(start_index>=0 && start_index < vps_poses_.size()){
            for(size_t i=start_index; i<vps_poses_.size(); ++i){
                if(vps_poses_.index(i).timestamp == timestamp){
                    printf("vps %d time %lf...\n", i, vps_poses_.index(i).timestamp);
                    return i;
                }
            }
        }
        return -1;
    }
    int VIOVPSFusion2::findTimeStampInVIO(int start_index, double timestamp){
        if(start_index>=0 && start_index < vio_poses_.size()){
            for(size_t i=start_index; i<vio_poses_.size(); ++i){
                if(vio_poses_.index(i).timestamp == timestamp){
                    printf("vio %d time %lf--%lf...\n", i, vio_poses_.index(i).timestamp, timestamp);
                    return i;
                }
            }
        }
        return -1;
    }

    void VIOVPSFusion2::MatchesFirstD1(){
        mlocker_.lock();
        if (viovps_matches_.index(0).first == 0)
        {
            printf("D1 pop front in viovps_matches_.\n");
            bool succ = viovps_matches_.pop_front();
            if(!succ){
                printf("no possible matches is zero...\n");
                exit(-1);
                return;
            }
        }
        for (int i = 0; i < viovps_matches_.size(); ++i)
        {
            viovps_matches_.index(i).first -= 1;
            if(viovps_matches_.index(i).first < 0){
                printf("d1:%d--", viovps_matches_.index(i).first);
                exit(-1);
            }
        }
        mlocker_.unlock();
    }
    void VIOVPSFusion2::MatchesSecondD1(){
        mlocker_.lock();
        if(viovps_matches_.index(0).second == 0){
            printf("D2 pop front in viovps_matches_.\n");
            bool succ = viovps_matches_.pop_front();
            if(!succ){
                printf("no possible matches is zero...\n");
                exit(-1);
                return;
            }
        }
        for(int i=0; i<viovps_matches_.size(); ++i){
            
            viovps_matches_.index(i).second -= 1;
            if(viovps_matches_.index(i).second < 0){
                printf("d2: %d\n", viovps_matches_.index(i).second);
                exit(-1);
            }
        }
        mlocker_.unlock();
    }
}   
}