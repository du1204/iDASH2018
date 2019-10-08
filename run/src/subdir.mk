################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/EncGWAS.cpp \
../src/Extractor.cpp \
../src/IDASH2018.cpp \
../src/TestGWAS.cpp

OBJS += \
./src/EncGWAS.o \
./src/Extractor.o \
./src/IDASH2018.o \
./src/TestGWAS.o

CPP_DEPS += \
./src/EncGWAS.d \
./src/Extractor.d \
./src/IDASH2018.d \
./src/TestGWAS.d

# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I/usr/local/include -I/home/kimandrik/HEAAN/HEAAN/src -O3 -g3 -Wall -c -fmessage-length=0 -std=c++11 -pthread -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


