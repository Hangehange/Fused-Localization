

#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "MyCorVTools.h"

#define MAXLINE 2501        //��ȡ����������

typedef struct {
    int id;                 //����
    double x, y, z;         //λ��xyz
    double utmX, utmY, utmZ; //���ϳ�ֵ����UTM����
    double qx, qy, qz, qw;  //��Ԫ��
    double rx, ry, rz, rw;  //�в�
}Lidar_pose_item;//�״����ݶ���

typedef struct {
    int id;                     //����
    double timestmp;            //ʱ���
    double utm_x, utm_y, utm_z; //utm����xyz
    double qx, qy, qz, qw;      //��Ԫ��
}INS_pose_item;//ins���ݶ���

typedef struct {
    double timestmp;            //ʱ���
    int sequence_num;           //���к�
    std::string sol_type;       //��λ״̬��ֻ��NARROW_INT��������
    double lat, lon, hig;       //γ�ȣ����ȣ��߶Ȳ���ֵ
    double utm_x, utm_y, utm_z; //ת���������
    double undulation;          //�������ֵ��ʵ�ʸ߶�;
    long double lat_dev, lon_dev, hig_dev;//ά�Ⱦ��ȸ߶ȵı�׼��

}GNSS_pose_item;//GNSS���ݶ���

typedef struct {
    double timestmp;                     //ʱ���
    long double lineA_x, lineA_y, lineA_z; //��������ϵxyz���ٶ���Ϣ
    long double AngV_x, AngV_y, AngV_z;  //��������ϵxyz���ٶ���Ϣ
}IMU_output_item;//IMU���ݽڵ㶨��

/*---------------------------------------------------------------------------------*/
/*
    ���ļ��ж�ȡ������ĺ���
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

    /*-----------------------дUTM��������----------------------*/

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
    std::string t;  //������������
    while (i <= _dataNum && fin.getline(line, sizeof(line))) {
        //fin.getline(line, sizeof(line));//��Ҫhead{
        fin.getline(line, sizeof(line));//Ҫ timestamp_sec: 1574149170.88
        std::stringstream timestamp_sec_word(line);     //����timestamp_sec��һ��
        timestamp_sec_word >> t;                        //  ����һ����
        timestamp_sec_word >> gnssList[i].timestmp; // ��ȡ��1574149170.88!
        fin.getline(line, sizeof(line));//��Ҫmodule_name: "gnss"
        fin.getline(line, sizeof(line));//��Ҫsequence_num: 23679
        fin.getline(line, sizeof(line));//��Ҫ}
        fin.getline(line, sizeof(line));//��Ҫmeasurement_time: 1258184388.8 
        fin.getline(line, sizeof(line));//��Ҫsol_status: SOL_COMPUTED
        fin.getline(line, sizeof(line));//Ҫsol_type: NARROW_INT 
        std::stringstream sol_type_word(line);          //����sol_type: NARROW_INT ��һ��
        sol_type_word >> t;                             //  ����һ����
        sol_type_word >> gnssList[i].sol_type;      //��ȡ��NARROW_INT ��
        fin.getline(line, sizeof(line));//Ҫlatitude: 28.1245135469 
        std::stringstream lat_word(line);              //����latitude: 28.1245135469 ��һ��
        lat_word >> t;                                 //  ����һ����
        lat_word >> gnssList[i].lat;               //��ȡ��28.1245135469��
        fin.getline(line, sizeof(line));//Ҫlongitude: 112.880556498 	
        std::stringstream lon_word(line);
        lon_word >> t;
        lon_word >> gnssList[i].lon;               //��ȡ��112.880556498 
        fin.getline(line, sizeof(line));//Ҫheight_msl: 58.0236365357 
        std::stringstream hig_word(line);
        hig_word >> t;
        hig_word >> gnssList[i].hig;               //��ȡ��58.0236365357 
        fin.getline(line, sizeof(line));//Ҫundulation: -16.8999996185 
        std::stringstream undul_word(line);
        undul_word >> t;
        undul_word >> gnssList[i].undulation;      //��ȡ��-16.8999996185
        fin.getline(line, sizeof(line));//��Ҫdatum_id: WGS84 
        fin.getline(line, sizeof(line));//Ҫlatitude_std_dev: 0.00532555021346
        std::stringstream lat_dev_word(line);
        lat_dev_word >> t;
        lat_dev_word >> gnssList[i].lat_dev;       //��ȡ��0.00532555021346
        fin.getline(line, sizeof(line));//Ҫlongitude_std_dev: 0.00600179005414
        std::stringstream lon_dev_word(line);
        lon_dev_word >> t;
        lon_dev_word >> gnssList[i].lon_dev;       //��ȡ��0.00600179005414
        fin.getline(line, sizeof(line));//Ҫheight_std_dev: 0.0142576200888 
        std::stringstream hi_dev_word(line);
        hi_dev_word >> t;
        hi_dev_word >> gnssList[i].hig_dev;       //��ȡ��0.0142576200888
        fin.getline(line, sizeof(line));//��Ҫbase_station_id: "0"
        fin.getline(line, sizeof(line));//��Ҫdifferential_age: 0.800000011921
        fin.getline(line, sizeof(line));//��Ҫsolution_age: 0.0
        fin.getline(line, sizeof(line));//��Ҫnum_sats_tracked: 30
        fin.getline(line, sizeof(line));//��Ҫnum_sats_in_solution: 30
        fin.getline(line, sizeof(line));//��Ҫnum_sats_l1: 30
        fin.getline(line, sizeof(line));//��Ҫnum_sats_multi: 26
        fin.getline(line, sizeof(line));//��Ҫextended_solution_status: 33
        fin.getline(line, sizeof(line));//��Ҫgalileo_beidou_used_mask: 48
        fin.getline(line, sizeof(line));//��Ҫgps_glonass_used_mask: 51
        fin.getline(line, sizeof(line));//��Ҫ����
        fin.getline(line, sizeof(line));//��Ҫ---
        i++;
    }
    fin.clear();
    fin.close();

    /*-------------------------------------WGS��UTM����---------------------------------------------*/

    WGS84Corr wgsTemp[MAXLINE];
    UTMCoor   utmTemp[MAXLINE];
    for (int i = 0; i < _dataNum; i++) {
        //��GNSSPostList����MAXLINE��WGS84Corr��UTMCoor
        wgsTemp[i].lat = gnssList[i].lat;
        wgsTemp[i].log = gnssList[i].lon;
        //��WGS84Corr����תUTMCoor
        utmTemp[i] = WGS_TO_UTM(wgsTemp[i]);
        gnssList[i].utm_x = utmTemp[i].x;
        gnssList[i].utm_y = utmTemp[i].y;
        gnssList[i].utm_z = gnssList[i].hig + gnssList[i].undulation;
        /*printf("%4d|%lf|%lf\n", i + 1, gnssList[i].utm_x, gnssList[i].utm_y);*/
    }
    //��MAXLINE��UTMCoorװ��GNSSPostList���utmȥ
}

void readTxt_IMU(IMU_output_item _imuList[],int _dataNum) {
    std::ifstream fin("imu_part.txt", std::ios::in);
    char line[1024] = { 0 };
    int i = 0;
    std::string t;  //������������
    while (i<= _dataNum && fin.getline(line, sizeof(line)) ) {
        fin.getline(line, sizeof(line));//Ҫ timestamp_sec: 1574149170.88
        std::stringstream timestamp_sec_word(line);     //����timestamp_sec��һ��
        timestamp_sec_word >> t;                        //  ����һ����
        timestamp_sec_word >> _imuList[i].timestmp; // ��ȡ��1574149170.88!
        fin.getline(line, sizeof(line));//��Ҫmodule_name: "gnss"
        fin.getline(line, sizeof(line));//��Ҫsequence_num: 23679
        fin.getline(line, sizeof(line));//��Ҫ}
        fin.getline(line, sizeof(line));//��Ҫmeasurement_time: 1258184388.89
        fin.getline(line, sizeof(line));//��Ҫmeasurement_span: 0.00800000037998
        fin.getline(line, sizeof(line));//��Ҫlinear_acceleration {
        fin.getline(line, sizeof(line));//Ҫx: 0.156353248026
        std::stringstream lax_word(line);
        lax_word >> t;
        lax_word >> _imuList[i].lineA_x;
        fin.getline(line, sizeof(line));//Ҫy: -0.0804513568985
        std::stringstream lay_word(line);
        lay_word >> t;
        lay_word >> _imuList[i].lineA_y;
        fin.getline(line, sizeof(line));//Ҫz: 9.42068214223
        std::stringstream laz_word(line);
        laz_word >> t;
        laz_word >> _imuList[i].lineA_z;
        fin.getline(line, sizeof(line));//��Ҫ}
        fin.getline(line, sizeof(line));//��Ҫangular_velocity {

        fin.getline(line, sizeof(line));//Ҫx: 0.156353248026
        std::stringstream anx_word(line);
        anx_word >> t;
        anx_word >> _imuList[i].AngV_x;
        fin.getline(line, sizeof(line));//Ҫy: -0.0804513568985
        std::stringstream any_word(line);
        any_word >> t;
        any_word >> _imuList[i].AngV_y;
        fin.getline(line, sizeof(line));//Ҫz: 9.42068214223
        std::stringstream anz_word(line);
        anz_word >> t;
        anz_word >> _imuList[i].AngV_z;
        fin.getline(line, sizeof(line));//��Ҫ}
        fin.getline(line, sizeof(line));//��Ҫ����
        fin.getline(line, sizeof(line));//��Ҫ��������
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