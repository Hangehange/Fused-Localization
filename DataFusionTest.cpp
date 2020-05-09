#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "txtLoadTool.h"

typedef struct {
    double fu_x, fu_y, fu_z;
    double fu_qx, fu_qy, fu_qz, fu_qw;
}Fusion_Pose_item;

/*---------------------------------------------------------------------------------*/
/*数据记录表*/
Lidar_pose_item  LidarPoseList[MAXLINE];

INS_pose_item    INSPoseList[MAXLINE];

GNSS_pose_item   GNSSPostList[MAXLINE];

IMU_output_item  IMUOutputList[MAXLINE];

Fusion_Pose_item FusionOutputList[MAXLINE]; //最终的融合输出


/*==================================核心：数据融合的方法===================================================*/
//计算Lidar浮动权重用的几个倍数参数
#define PARA_QN  3
#define PARA_XYZ 4
/*宏定义的这几个数字都要去datasheet的histogram里看着图调*/

//四元数残差的阈值宏定义
#define THRESHOLD_QX 0.0652
#define THRESHOLD_QY 0.0650
#define THRESHOLD_QZ 0.07
#define THRESHOLD_QW 0.046
//Lidar与INS的差阈值定义
#define THRESHOLD_DELTA_QX 0.032
#define THRESHOLD_DELTA_QY 0.0125
#define THRESHOLD_DELTA_QZ 1.98
#define THRESHOLD_DELTA_QW 1.8
#define THRESHOLD_DELTA_X 16.5
#define THRESHOLD_DELTA_Y 17.5
#define THRESHOLD_DELTA_Z 1.7
//滤波器窗口大小
#define WINDOW_WIDTH_MID 5
#define WINDOW_WIDTH_AVG 3

//定义计算步骤（计算四元数还是计算三维坐标）
#define STEP_QN     1234
#define STEP_XYZ    4321
//传感器状态异常失灵错误码
#define STATE_ERROR 0

typedef struct {
    double qx_weight, qy_weight, qz_weight, qw_weight;
}LidarQNweight;

typedef struct{
    double x_weight,y_weight,z_weight;
}GnssXYZWeight;

typedef struct {
    double INS_Bweight, GNSS_Bweight, Lidar_Bweight;//基础权重
    LidarQNweight lidarQnWeight;//加上浮动权重之后的LIDAR总权重
    GnssXYZWeight gnssXYZweight;//加上浮动权重之后的GNSS总权重
}WeightStructDef;//定义一则数据的各个权重

WeightStructDef weightStruct ;//权重的记录表

typedef struct {
    int LidarState;//整体状态：若四个里有两个超过阈值就是失灵
    int qxState, qyState, qzState, qwState;//正常为1，失灵为0
}LidarStateDef; //雷达状态的结构体

typedef int GNSSStateDef;             //GNSS状态定义（实际就是个int）

/*   √定量GNSS精度。 NARROW_INT=1；NARROW_FLOAT=2;PSRDIFF=0;*/
GNSSStateDef getGNSStype(GNSS_pose_item _item) {
    if (_item.sol_type == "NARROW_INT")        return 1;
    else if (_item.sol_type == "PSRDIFF" 
        or _item.sol_type=="SINGLE")           return STATE_ERROR;
    else if (_item.sol_type == "NARROW_FLOAT") return 2;
    else return -1;//错误码
}

/*   √输入一条Lidar记录返回一个LidarStateDef结构体,其中包括各个值是否超过阈值的判断 */
LidarStateDef getLidarState(Lidar_pose_item _item) {
    LidarStateDef stateTemp = { 1,1,1,1,1 };//初始状态默认都正常
    int count = 0;
    if (_item.rx > THRESHOLD_QX) {
        stateTemp.qxState = STATE_ERROR;
        count++;
    }
    if (_item.ry > THRESHOLD_QY) {
        stateTemp.qyState = STATE_ERROR;
        count++;
    }
    if (_item.rz > THRESHOLD_QZ) {
        stateTemp.qzState = STATE_ERROR;
        count++;
    }
    if (_item.rw > THRESHOLD_QW) {
        stateTemp.qwState = STATE_ERROR;
        count++;
    }
    if (count >= 2)
        stateTemp.LidarState = STATE_ERROR;
    return stateTemp;
}

/*   √Lidar和INS相互校验。计算INS和Lidar测量四元数的差，并和差阈值相比较
        inputs:分别传进来Lidar和INS的一条记录
        return:计算的置信度conf
*/
double insLidarMutualCheck(Lidar_pose_item _lidarItem, INS_pose_item _insItem) {
    double conf;//整体置信度
    double confX, confY, confZ, confQX, confQY, confQZ,confQW;
    confQX = fabs(fabs(_lidarItem.qx - _insItem.qx) - THRESHOLD_DELTA_QX) / THRESHOLD_DELTA_QX;
    confQY = fabs(fabs(_lidarItem.qy - _insItem.qy) - THRESHOLD_DELTA_QY) / THRESHOLD_DELTA_QY;
    confQZ = fabs(fabs(_lidarItem.qz - _insItem.qz) - THRESHOLD_DELTA_QZ) / THRESHOLD_DELTA_QZ;
    confQW = fabs(fabs(_lidarItem.qw - _insItem.qw) - THRESHOLD_DELTA_QW) / THRESHOLD_DELTA_QW;
    confX = fabs(fabs(_lidarItem.utmX - _insItem.utm_x) - THRESHOLD_DELTA_X) / THRESHOLD_DELTA_X;
    confY = fabs(fabs(_lidarItem.utmY - _insItem.utm_y) - THRESHOLD_DELTA_Y) / THRESHOLD_DELTA_Y;
    confZ = fabs(fabs(_lidarItem.utmZ - _insItem.utm_z) - THRESHOLD_DELTA_Z) / THRESHOLD_DELTA_Z;
    conf = (confX + confY + confZ + confQX + confQY + confQZ + confQW) / 7;
    return conf;
}

/*    √核心方法：按照不同的精度要求给预设权重赋值
        para:分别传进来Lidar和Gnss的一条记录进来，
        step表示的是计算四元数的步骤还是计算坐标步骤: STEP_QN 或者STEP_XYZ
*/
void basicWeightLoad(Lidar_pose_item _lidarItem, GNSS_pose_item _gnssItem,INS_pose_item _insItem,int _step) {
    weightStruct.INS_Bweight = weightStruct.GNSS_Bweight = weightStruct.Lidar_Bweight = 5;//先默认基础权重都给5
    GNSSStateDef  Gstate = getGNSStype(_gnssItem);
    LidarStateDef Lstate = getLidarState(_lidarItem);
    if      (STEP_QN == _step) //计算四元数环节给权重赋值
    {    
        if (STATE_ERROR == Lstate.LidarState )//先Lidar状态正确:检查_lidarItem的LidarState
        {
            //为0说明Lidar失灵，Lidar至少有两个残差超阈值
            weightStruct.Lidar_Bweight -= 2;
            weightStruct.INS_Bweight += 1;
        }
        else 
        {
            switch (getGNSStype(_gnssItem)) {//先看GNSS状态
                case 1: //GNSS状态NARROW_INT，基础权值不动，INS没有浮动权重项，
                {
                    break;
                }
                case 2: {//GNSS状态NARROW_FLOAT
                    weightStruct.Lidar_Bweight += 1;
                    weightStruct.INS_Bweight -= 1;//按表规定调节基础权值

                    break;
                }
                case 0: {//GNSS状态PSRIFF
                    weightStruct.Lidar_Bweight += 3;
                    weightStruct.INS_Bweight -= -1;//按表规定调节基础权值
                    //GNSS状态不好，INS可信度下降，INS总权重 = INS预设权重 * conf，
                    weightStruct.INS_Bweight *= insLidarMutualCheck(_lidarItem, _insItem);
                    break;
                }
                case -1: {//错误
                    printf("error in GNSS type\n");
                    break;
                }
                default: {
                    break;
                }
            }
        }     
        // 调不用分类的Lidar浮动权重：Lidar权重Ws[4] = Lidar浮动权重w[4] - 残差[4] * PARA
        weightStruct.lidarQnWeight.qx_weight =  weightStruct.Lidar_Bweight - _lidarItem.rx * PARA_QN;
        if (weightStruct.lidarQnWeight.qx_weight < 0)
            weightStruct.lidarQnWeight.qx_weight = 0;
        weightStruct.lidarQnWeight.qy_weight =  weightStruct.Lidar_Bweight - _lidarItem.ry * PARA_QN;
        if (weightStruct.lidarQnWeight.qy_weight < 0)
            weightStruct.lidarQnWeight.qy_weight = 0;
        weightStruct.lidarQnWeight.qz_weight =  weightStruct.Lidar_Bweight - _lidarItem.rz * PARA_QN;
        if (weightStruct.lidarQnWeight.qz_weight < 0)
            weightStruct.lidarQnWeight.qz_weight = 0;
        weightStruct.lidarQnWeight.qw_weight =  weightStruct.Lidar_Bweight - _lidarItem.rw * PARA_QN;
        if (weightStruct.lidarQnWeight.qw_weight < 0)
            weightStruct.lidarQnWeight.qw_weight = 0;
    }
    //计算三维坐标环节给权重赋值
    else if (STEP_XYZ == _step) //计算坐标环节的给权重赋值
    {
        if (STATE_ERROR == Lstate.LidarState)//先Lidar状态正确:检查_lidarItem的LidarState
        {
            //为0说明Lidar失灵，Lidar至少有两个残差超阈值
            weightStruct.Lidar_Bweight -= 2;
            weightStruct.INS_Bweight   += 1;
            weightStruct.GNSS_Bweight  += 1;
        }
        else //Lidar正常
        {          
            switch (getGNSStype(_gnssItem)) {//先看GNSS状态
                case 1: //INS没有浮动权重项，预设权重就是总权重， Lidar浮动权重w[4] = conf * Lidar固定权重
                {
                    double confTemp = insLidarMutualCheck(_lidarItem, _insItem); 
                    weightStruct.lidarQnWeight.qx_weight 
                        = weightStruct.Lidar_Bweight - _lidarItem.rx * PARA_QN * confTemp;
                    weightStruct.lidarQnWeight.qy_weight 
                        = weightStruct.Lidar_Bweight - _lidarItem.ry * PARA_QN * confTemp;
                    weightStruct.lidarQnWeight.qz_weight 
                        = weightStruct.Lidar_Bweight - _lidarItem.rz * PARA_QN * confTemp;
                    weightStruct.lidarQnWeight.qw_weight 
                        = weightStruct.Lidar_Bweight - _lidarItem.rw * PARA_QN * confTemp;
                    break;
                }
                case 2: {//GNSS状态NARROW_FLOAT
                    weightStruct.Lidar_Bweight  += 2;
                    weightStruct.INS_Bweight    -= 1;
                    weightStruct.GNSS_Bweight   -= 2;//按表规定调节基础权值
                    break;
                }
                case 0: {//GNSS状态PSRIFF，INS总权重 = INS预设权重 * conf， Lidar浮动权重w[4] = para * 残差
                    weightStruct.Lidar_Bweight += 3;
                    weightStruct.GNSS_Bweight  -= 3;  
                    weightStruct.INS_Bweight   -= 1; //按表规定调节基础权值
                    //GNSS状态不好，INS可信度下降，INS总权重 = INS预设权重 * conf，
                    weightStruct.INS_Bweight *= insLidarMutualCheck(_lidarItem, _insItem);
                    break;
                }
                case -1: {//错误
                    printf("error in GNSS type\n");
                    break;
                }
                default: {
                    break;
                }
            }
            weightStruct.Lidar_Bweight *= insLidarMutualCheck(_lidarItem, _insItem);
            //最后设定GNSS权重Ws[3] = W - std_dev[3] * para 
            //若出现负值，说明差的太远了，直接扔掉
            weightStruct.gnssXYZweight.x_weight = weightStruct.GNSS_Bweight - _gnssItem.lat_dev * PARA_XYZ;
            if (weightStruct.gnssXYZweight.x_weight < 0)
                weightStruct.gnssXYZweight.x_weight = 0;
            weightStruct.gnssXYZweight.y_weight = weightStruct.GNSS_Bweight - _gnssItem.lon_dev * PARA_XYZ;
            if (weightStruct.gnssXYZweight.y_weight < 0)
                weightStruct.gnssXYZweight.y_weight = 0;
            weightStruct.gnssXYZweight.z_weight = weightStruct.GNSS_Bweight - _gnssItem.hig_dev * PARA_XYZ;
            if (weightStruct.gnssXYZweight.z_weight < 0)
                weightStruct.gnssXYZweight.z_weight = 0;
        }
    }
    else {
        printf("step code wrong!");
    }
}



/*滑动平均实现平滑处理:
    将滑动窗口大小windowSize内的数据的算术平均值作为_curr的值
*/
Fusion_Pose_item averageFiltering(Fusion_Pose_item _fuList[],int _currID)
{
    int start = 0;//窗口的开始和结束
    if (_currID == 0) {
        return _fuList[0];
    }
    else if (_currID >= WINDOW_WIDTH_AVG)
    {
        start = _currID - WINDOW_WIDTH_AVG;
    }
    double sumX, sumY, sumZ, sumQX, sumQY, sumQZ, sumQW;
    sumQW = sumQZ = sumQY = sumQX = sumX = sumY = sumZ = 0;
    for (int i = _currID; i > start; i--)
    {
        sumX += _fuList[i].fu_x;     sumY += _fuList[i].fu_y;  sumZ += _fuList[i].fu_z;
        sumQX += _fuList[i].fu_qx;   sumQY += _fuList[i].fu_qy;
        sumQZ += _fuList[i].fu_qz;   sumQW += _fuList[i].fu_qw;
    }
    Fusion_Pose_item result;
    result.fu_x = sumX / ((double)_currID - (double)start);
    result.fu_y = sumY / ((double)_currID - (double)start);
    result.fu_z = sumZ / ((double)_currID - (double)start);
    result.fu_qx = sumQX / ((double)_currID - (double)start);
    result.fu_qy = sumQY / ((double)_currID - (double)start);
    result.fu_qz = sumQZ / ((double)_currID - (double)start);
    result.fu_qw = sumQW / ((double)_currID - (double)start);

    return result;
    

}

/*简单选择排序*/
#define SWAP(x, y, t)  ((t) = (x), (x) = (y), (y) = (t))
void sortSquence(double list[], int n)
{
    int i, j, min;
    double temp;
    for (i = 0; i < n - 1; i++) {
        min = i;
        for (j = i + 1; j < n; j++)
            if (list[j] < list[min])
                min = j;
        SWAP(list[i], list[min], temp);
    }
}


/* 中值滤波实现平滑处理
    将奇数滑动窗口大小windowSize内的数据排序，以中位数当作输出
    step:   1.currentID左边右边各自为0.5*(_windowSize-1)，以此作为window
            2.选择排序
    */

Fusion_Pose_item medianFilteringg(Fusion_Pose_item _fuList[], int _currID)
{
    Fusion_Pose_item result;
    int halfWindow = (WINDOW_WIDTH_MID - 1)/2; //半窗口长度
    double tempX[WINDOW_WIDTH_MID], tempY[WINDOW_WIDTH_MID], tempZ[WINDOW_WIDTH_MID];
    double tempQX[WINDOW_WIDTH_MID], tempQY[WINDOW_WIDTH_MID], tempQZ[WINDOW_WIDTH_MID], tempQW[WINDOW_WIDTH_MID];
    if (_currID <= WINDOW_WIDTH_MID or _fuList[_currID+ halfWindow].fu_x == 0) //最开始的和最末尾的不处理
    {
        return _fuList[_currID];
    }
    else {
        for (int i = _currID - halfWindow; i < _currID + halfWindow; i++) {
            tempX[i] = _fuList[i].fu_x; tempY[i] = _fuList[i].fu_y; tempZ[i] = _fuList[i].fu_z;//给他们赋值
            tempQX[i] = _fuList[i].fu_qx; tempQY[i] = _fuList[i].fu_qy;
            tempQZ[i] = _fuList[i].fu_qz; tempQW[i] = _fuList[i].fu_qw;
        }
        sortSquence(tempX , WINDOW_WIDTH_MID); sortSquence(tempY , WINDOW_WIDTH_MID);  sortSquence(tempZ, WINDOW_WIDTH_MID);
        sortSquence(tempQX, WINDOW_WIDTH_MID); sortSquence(tempQY, WINDOW_WIDTH_MID);
        sortSquence(tempQZ, WINDOW_WIDTH_MID); sortSquence(tempQW, WINDOW_WIDTH_MID);
    }
    result = { tempX[_currID],tempY[_currID],tempZ[_currID],tempQX[_currID],tempQY[_currID],tempQZ[_currID],tempQW[_currID] };
    return result;
}


/*加权平均:*/
double weightAverage(double weight[], double data[],int _lenth)
{
    double sumWeight = 0; double sumResult = 0;
    for (size_t i = 0; i < _lenth; i++)
    {
        sumResult += data[i] * weight[i];
        sumWeight += weight[i];
    }
    return sumResult / sumWeight;
}

/*核心方法：一则数据数据相互检验相互融合*/
Fusion_Pose_item weightMean(Lidar_pose_item _Litem, GNSS_pose_item _Gitem, INS_pose_item _Iitem) {
    Fusion_Pose_item result;
    WeightStructDef tempW;
    //第一步：先加权平均Lidar和INS计算四元数
    tempW = weightStruct;    
    basicWeightLoad(_Litem, _Gitem, _Iitem, STEP_QN);//计算四元数的权值weightStruct
    result.fu_qx = (tempW.lidarQnWeight.qx_weight * _Litem.qx + tempW.INS_Bweight * _Iitem.qx)
        / (tempW.lidarQnWeight.qx_weight + tempW.INS_Bweight);
    result.fu_qy = tempW.lidarQnWeight.qy_weight * _Litem.qy + tempW.INS_Bweight * _Iitem.qy
        / (tempW.lidarQnWeight.qy_weight + tempW.INS_Bweight);
    result.fu_qz = tempW.lidarQnWeight.qz_weight * _Litem.qz + tempW.INS_Bweight * _Iitem.qz
        / (tempW.lidarQnWeight.qz_weight + tempW.INS_Bweight);
    result.fu_qw = tempW.lidarQnWeight.qw_weight * _Litem.qw + tempW.INS_Bweight * _Iitem.qw
        / (tempW.lidarQnWeight.qw_weight + tempW.INS_Bweight);
    //第二步：加权平均三个传感器计算XYZ坐标
    basicWeightLoad(_Litem, _Gitem, _Iitem, STEP_XYZ);//计算四元数的权值weightStruct
    tempW = weightStruct;
    result.fu_x = (tempW.gnssXYZweight.x_weight * _Gitem.utm_x + tempW.INS_Bweight * _Iitem.utm_x
        + tempW.Lidar_Bweight * _Litem.utmX) / (tempW.gnssXYZweight.x_weight + tempW.INS_Bweight + tempW.Lidar_Bweight);
    result.fu_y = (tempW.gnssXYZweight.y_weight * _Gitem.utm_y + tempW.INS_Bweight * _Iitem.utm_y
        + tempW.Lidar_Bweight * _Litem.utmY) / (tempW.gnssXYZweight.y_weight + tempW.INS_Bweight + tempW.Lidar_Bweight);
    result.fu_z = (tempW.gnssXYZweight.z_weight * _Gitem.utm_z + tempW.INS_Bweight * _Iitem.utm_z
        + tempW.Lidar_Bweight * _Litem.utmZ) / (tempW.gnssXYZweight.z_weight + tempW.INS_Bweight + tempW.Lidar_Bweight);
    return result;
}

/*按步骤融合数据*/
void finalDataOutput(int size) {
    //这一个循环相当于一辆车在行驶过程中，接受源源不断的流数据，将结果写道FusionOutputList[i]中
    for (int i = 0; i < size; i++) 
    {
        FusionOutputList[i] = weightMean(LidarPoseList[i], GNSSPostList[i], INSPoseList[i]);//计算加权平均
        FusionOutputList[i] = medianFilteringg(FusionOutputList, i);                        //中值滤波
        FusionOutputList[i] = averageFiltering(FusionOutputList, i);                        //均值滤波
    }
}

int main()
{
    int dataSize = 15;  //一次读取多少条数据
    if (dataSize > 2500) {
        printf("没那么多\n");
        return -1;
    }
    else {
        loadAllDataFromTxt(GNSSPostList, INSPoseList, LidarPoseList, IMUOutputList, dataSize);
        finalDataOutput(dataSize);
        for (int i = 1; i < dataSize; i++) {
            printf("%3d| %lf|%lf|%lf|%lf|%lf|%lf|%lf|\n", i,
                FusionOutputList[i].fu_x, FusionOutputList[i].fu_y, FusionOutputList[i].fu_z,
                FusionOutputList[i].fu_qx, FusionOutputList[i].fu_qy,
                FusionOutputList[i].fu_qw, FusionOutputList[i].fu_qz);
        }
        return 0;
    }
}