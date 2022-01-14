#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <string>
#include "src/utils/basic_types.h"
#include "src/core/svvfusion.h"
#include <unistd.h>
// model just support euroc
void readposecsv(std::string filepath, std::vector<svv_fusion::Posed_t>& poses){
    std::ifstream reader(filepath);
    poses.reserve(5000);
    if(reader.is_open()){
        std::string line;
        std::getline(reader, line);
        while (std::getline(reader, line))
        {
            // split by ","
            std::replace(line.begin(), line.end(), ',', ' ');
            std::stringstream ss(line);
            svv_fusion::Posed_t p;

            ss >> p.timestamp;
            ss >> p.t_wc[0];
            ss >> p.t_wc[1];
            ss >> p.t_wc[2];
            ss >> p.q_wc.w();
            ss >> p.q_wc.x();
            ss >> p.q_wc.y();
            ss >> p.q_wc.z();

            poses.push_back(p);
        }
        reader.close();
    }
    else{
        printf("cannot open the file %s\n", filepath.c_str());
        exit(-1);
    }
}

int main(int argc, char* argv[]){
    // 把 两个pose 数据读取出来，分别push into svvfusion
    std::vector<svv_fusion::Posed_t> vioposes, vpsposes;
    readposecsv("/home/lut/Desktop/evo/vio.csv", vioposes);
    readposecsv("/home/lut/Desktop/evo/vps.csv", vpsposes);

    std::vector<svv_fusion::Posed_t> fusion_results;
    svv_fusion::VIOVPSFusion vvfusion(50, 150, 50);
    int vioposition = 0, vpsposition = 0;
    while(true){
        if(vioposes.size() <= vioposition || vpsposes.size() <= vpsposition){
            break;
        }
        auto viot = vioposes[vioposition].timestamp;
        auto vpst = vpsposes[vpsposition].timestamp;
        // printf("the vio/vps pose: %lf, %lf...\n", viot, vpst);
        // std::cout<<"push in viop: "<<vioposes[vioposition]<<std::endl;
        // std::cout<<"push in vpsp: "<<vpsposes[vpsposition]<<std::endl;

        if(!vioposes[vioposition].valued() || !vpsposes[vpsposition].valued()){
            printf("data error...\n");
            exit(-1);
        }
        if(viot == vpst){
            printf("push vps vio.\n");
            vvfusion.PushVIOPose(vioposes[vioposition]);
            vvfusion.PushVPSPose(vpsposes[vpsposition]);           
            ++vioposition;
            ++vpsposition;
        }
        else if(viot < vpst){
            printf("push vio.\n");
            vvfusion.PushVIOPose(vioposes[vioposition]);
            ++vioposition;
        }
        else if(viot > vpst){
            printf("push vps vio.\n");
            vvfusion.PushVPSPose(vpsposes[vpsposition]);
            ++vpsposition;
        }
        usleep(10000);
    }
    return 0;
}