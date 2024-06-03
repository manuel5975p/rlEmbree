#include <embree4/rtcore.h>
#include <stdio.h>

int main(int argc, char* argv[]){
    RTCDevice dev;
    dev = rtcNewDevice(NULL);
    printf("Helo\n");
}