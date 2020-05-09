

#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "MyCorVTools.h"

#define MAXLINE 2501        //读取多少行数据

typedef struct {
    int id;                 //索引
    double x, y, z;         //位置xyz
    double utmX, utmY, utmZ; //加上初值的真UTM坐标
    double qx, qy, qz, qw;  //四元数
    double rx, ry, rz, rw;  //残差
}Lidar_pose_item;//雷达数据定义

typedef struct {
    int id;                     //索引
    double timestmp;            //时间戳
    double utm_x, utm_y, utm_z; //utm坐标xyz
    double qx, qy, qz, qw;      //四元数
}INS_pose_item;//ins数据定义

typedef struct {
    double timestmp;            //时间戳
    int sequence_num;           //序列号
    std::string sol_type;       //定位状态：只有NARROW_INT才是良好
    double lat, lon, hig;       //纬度，经度，高度测量值
    double utm_x, utm_y, utm_z; //转换后的坐标
    double undulation;          //加上这个值是实际高度;
    long double lat_dev, lon_dev, hig_dev;//维度经度高度的标准差

}GNSS_pose_item;//GNSS数据定义

typedef struct {
    double timestmp;                     //时间戳
    long double lineA_x, lineA_y, lineA_z; //车身坐标系xyz加速度信息
    long double AngV_x, AngV_y, AngV_z;  //车身坐标系xyz角速度信息
}IMU_output_item;//IMU数据节点定义

/*---------------------------------------------------------------------------------*/
/*
    从文件中读取到程序的函数
*/
void readTxt_INS(INS_pose_item insList[], int _dataNum) {
    std::ifstream fin("ins_pose.txt", std::ios::in);
    char line[1024] = { 0 };
    int i = 0;
    while (i <= _dataNum && fin.getline(line, sizeof(line)))
    {
        std::stringstream word(line);
        word >> insList[i].id; word >> insList[i].timestmp;
        word >> insList[i].utm_x; word >> insList[i].utm_y; word >> insList[i].utm_z;
        word >> insList[i].qx; word >> insList[i].qy; word >> insList[i].qz; word >> insList[i].qw;
        i++;
    }
    fin.clear();
    fin.close();
}

void readTxt_Lidar(Lidar_pose_item _lidarList[], int _dataNum) {
    std::ifstream fin("lidar_pose.txt", std::ios::in);
    char line[1024] = { 0 };
    int i = 0;
    while (i <= _dataNum && fin.getline(line, sizeof(line)))
    {
        std::stringstream word(line);
        word >> _lidarList[i].id;
        word >> _lidarList[i].qx; word >> _lidarList[i].qy; word >> _lidarList[i].qz; word >> _lidarList[i].qw;
        word >> _lidarList[i].x;  word >> _lidarList[i].y;  word >> _lidarList[i].z;
        word >> _lidarList[i].rx; word >> _lidarList[i].ry; word >> _lidarList[i].rz; word >> _lidarList[i].rw;
        i++;
    }
    fin.clear();
    fin.close();

    /*-----------------------写UTM三个坐标----------------------*/

    for (i = 0; i < _dataNum; i++) {
        _lidarList[i].utmX = _lidarList[i].x + 684704.4281840;
        _lidarList[i].utmY = _lidarList[i].y + 3112421.914861;
        _lidarList[i].utmZ = _lidarList[i].z + 41.810467;

    }
}

void readTxt_GNSS(GNSS_pose_item gnssList[],int _dataNum) {
    std::ifstream fin("gnss_best_pose.txt", std::ios::in);
    char line[1024] = { 0 };
    int i = 0;
    std::string t;  //用来倒垃圾的
    while (i <= _dataNum && fin.getline(line, sizeof(line))) {
        //fin.getline(line, sizeof(line));//不要head{
        fin.getline(line, sizeof(line));//要 timestamp_sec: 1574149170.88
        std::stringstream timestamp_sec_word(line);     //进入timestamp_sec这一行
        timestamp_sec_word >> t;                        //  往右一个词
        timestamp_sec_word >> gnssList[i].timestmp; // 读取到1574149170.88!
        fin.getline(line, sizeof(line));//不要module_name: "gnss"
        fin.getline(line, sizeof(line));//不要sequence_num: 23679
        fin.getline(line, sizeof(line));//不要}
        fin.getline(line, sizeof(line));//不要measurement_time: 1258184388.8 
        fin.getline(line, sizeof(line));//不要sol_status: SOL_COMPUTED
        fin.getline(line, sizeof(line));//要sol_type: NARROW_INT 
        std::stringstream sol_type_word(line);          //进入sol_type: NARROW_INT 这一行
        sol_type_word >> t;                             //  往右一个词
        sol_type_word >> gnssList[i].sol_type;      //读取到NARROW_INT ！
        fin.getline(line, sizeof(line));//要latitude: 28.1245135469 
        std::stringstream lat_word(line);              //进入latitude: 28.1245135469 这一行
        lat_word >> t;                                 //  往右一个词
        lat_word >> gnssList[i].lat;               //读取到28.1245135469！
        fin.getline(line, sizeof(line));//要longitude: 112.880556498 	
        std::stringstream lon_word(line);
        lon_word >> t;
        lon_word >> gnssList[i].lon;               //读取到112.880556498 
        fin.getline(line, sizeof(line));//要height_msl: 58.0236365357 
        std::stringstream hig_word(line);
        hig_word >> t;
        hig_word >> gnssList[i].hig;               //读取到58.0236365357 
        fin.getline(line, sizeof(line));//要undulation: -16.8999996185 
        std::stringstream undul_word(line);
        undul_word >> t;
        undul_word >> gnssList[i].undulation;      //读取到-16.8999996185
        fin.getline(line, sizeof(line));//不要datum_id: WGS84 
        fin.getline(line, sizeof(line));//要latitude_std_dev: 0.00532555021346
        std::stringstream lat_dev_word(line);
        lat_dev_word >> t;
        lat_dev_word >> gnssList[i].lat_dev;       //读取到0.00532555021346
        fin.getline(line, sizeof(line));//要longitude_std_dev: 0.00600179005414
        std::stringstream lon_dev_word(line);
        lon_dev_word >> t;
        lon_dev_word >> gnssList[i].lon_dev;       //读取到0.00600179005414
        fin.getline(line, sizeof(line));//要height_std_dev: 0.0142576200888 
        std::stringstream hi_dev_word(line);
        hi_dev_word >> t;
        hi_dev_word >> gnssList[i].hig_dev;       //读取到0.0142576200888
        fin.getline(line, sizeof(line));//不要base_station_id: "0"
        fin.getline(line, sizeof(line));//不要differential_age: 0.800000011921
        fin.getline(line, sizeof(line));//不要solution_age: 0.0
        fin.getline(line, sizeof(line));//不要num_sats_tracked: 30
        fin.getline(line, sizeof(line));//不要num_sats_in_solution: 30
        fin.getline(line, sizeof(line));//不要num_sats_l1: 30
        fin.getline(line, sizeof(line));//不要num_sats_multi: 26
        fin.getline(line, sizeof(line));//不要extended_solution_status: 33
        fin.getline(line, sizeof(line));//不要galileo_beidou_used_mask: 48
        fin.getline(line, sizeof(line));//不要gps_glonass_used_mask: 51
        fin.getline(line, sizeof(line));//不要空行
        fin.getline(line, sizeof(line));//不要---
        i++;
    }
    fin.clear();
    fin.close();

    /*-------------------------------------WGS换UTM部分---------------------------------------------*/

    WGS84Corr wgsTemp[MAXLINE];
    UTMCoor   utmTemp[MAXLINE];
    for (int i = 0; i < _dataNum; i++) {
        //拿GNSSPostList构造MAXLINE个WGS84Corr和UTMCoor
        wgsTemp[i].lat = gnssList[i].lat;
        wgsTemp[i].log = gnssList[i].lon;
        //把WGS84Corr挨个转UTMCoor
        utmTemp[i] = WGS_TO_UTM(wgsTemp[i]);
        gnssList[i].utm_x = utmTemp[i].x;
        gnssList[i].utm_y = utmTemp[i].y;
        gnssList[i].utm_z = gnssList[i].hig + gnssList[i].undulation;
        /*printf("%4d|%lf|%lf\n", i + 1, gnssList[i].utm_x, gnssList[i].utm_y);*/
    }
    //把MAXLINE个UTMCoor装到GNSSPostList里的utm去
}

void readTxt_IMU(IMU_output_item _imuList[],int _dataNum) {
    std::ifstream fin("imu_part.txt", std::ios::in);
    char line[1024] = { 0 };
    int i = 0;
    std::string t;  //用来倒垃圾的
    while (i<= _dataNum && fin.getline(line, sizeof(line)) ) {
        fin.getline(line, sizeof(line));//要 timestamp_sec: 1574149170.88
        std::stringstream timestamp_sec_word(line);     //进入timestamp_sec这一行
        timestamp_sec_word >> t;                        //  往右一个词
        timestamp_sec_word >> _imuList[i].timestmp; // 读取到1574149170.88!
        fin.getline(line, sizeof(line));//不要module_name: "gnss"
        fin.getline(line, sizeof(line));//不要sequence_num: 23679
        fin.getline(line, sizeof(line));//不要}
        fin.getline(line, sizeof(line));//不要measurement_time: 1258184388.89
        fin.getline(line, sizeof(line));//不要measurement_span: 0.00800000037998
        fin.getline(line, sizeof(line));//不要linear_acceleration {
        fin.getline(line, sizeof(line));//要x: 0.156353248026
        std::stringstream lax_word(line);
        lax_word >> t;
        lax_word >> _imuList[i].lineA_x;
        fin.getline(line, sizeof(line));//要y: -0.0804513568985
        std::stringstream lay_word(line);
        lay_word >> t;
        lay_word >> _imuList[i].lineA_y;
        fin.getline(line, sizeof(line));//要z: 9.42068214223
        std::stringstream laz_word(line);
        laz_word >> t;
        laz_word >> _imuList[i].lineA_z;
        fin.getline(line, sizeof(line));//不要}
        fin.getline(line, sizeof(line));//不要angular_velocity {

        fin.getline(line, sizeof(line));//要x: 0.156353248026
        std::stringstream anx_word(line);
        anx_word >> t;
        anx_word >> _imuList[i].AngV_x;
        fin.getline(line, sizeof(line));//要y: -0.0804513568985
        std::stringstream any_word(line);
        any_word >> t;
        any_word >> _imuList[i].AngV_y;
        fin.getline(line, sizeof(line));//要z: 9.42068214223
        std::stringstream anz_word(line);
        anz_word >> t;
        anz_word >> _imuList[i].AngV_z;
        fin.getline(line, sizeof(line));//不要}
        fin.getline(line, sizeof(line));//不要换行
        fin.getline(line, sizeof(line));//不要————
        i++;
    }
    //printf("%lf", IMUOutputList[i].AngV_z);
    fin.clear();
    fin.close();
}

void loadAllDataFromTxt(GNSS_pose_item gnssList[], INS_pose_item insList[], Lidar_pose_item lidarList[],IMU_output_item imuList[],int _dataNum) 
{

    readTxt_IMU(imuList, _dataNum);  
    readTxt_GNSS(gnssList, _dataNum);  
    readTxt_INS(insList, _dataNum);    
    readTxt_Lidar(lidarList, _dataNum);  
}