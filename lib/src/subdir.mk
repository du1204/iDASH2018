################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/CyclicArray.cpp \
../src/EncGWAS.cpp \
../src/Extractor.cpp \
../src/IDASH2018.cpp \
../src/PlainGWAS.cpp \
../src/SchemeGD.cpp \
../src/TestCyclicArray.cpp \
../src/TestGD.cpp \
../src/TestGWAS.cpp \
../src/TrueGWAS.cpp \

OBJS += \
./src/CyclicArray.o \
./src/EncGWAS.o \
./src/Extractor.o \
./src/IDASH2018.o \
./src/PlainGWAS.o \
./src/SchemeGD.o \
./src/TestCyclicArray.o \
./src/TestGD.o \
./src/TestGWAS.o \
./src/TrueGWAS.o \


CPP_DEPS += \
./src/CyclicArray.d \
./src/EncGWAS.d \
./src/Extractor.d \
./src/IDASH2018.d \
./src/PlainGWAS.d \
./src/SchemeGD.d \
./src/TestCyclicArray.d \
./src/TestGD.d \
./src/TestGWAS.d \
./src/TrueGWAS.d \


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/Users/KimDuhyeong/Public/Lib/FASTHEAAN/FASTHEAAN/src -O3 -pthread -c -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


