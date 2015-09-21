compile:
	@g++ -c src_mn/MnTagClust.cpp src_mn/MnNucPos.cpp src_mn/MnMain.cpp src_ho/HoTagClust.cpp src_ho/HoNucPos.cpp src_ho/HoMain.cpp src_common/Main.cpp
	@g++ -o numster MnTagClust.o MnNucPos.o MnMain.o HoTagClust.o HoNucPos.o HoMain.o Main.o -pthread
	@rm *.o

