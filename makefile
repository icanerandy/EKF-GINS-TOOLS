gen:
	mkdir -p build/
	conan install . -s build_type=Debug --build=missing
	conan install . -s build_type=Release --build=missing
	conan install . -s build_type=RelWithDebInfo --build=missing
	conan install . -s build_type=MinSizeRel --build=missing
	cmake --preset conan-default --log-level=VERBOSE

md:
	cmake --build --preset conan-debug --parallel

mrd:
	cmake --build --preset conan-relwithdebinfo --parallel

mr:
	cmake --build --preset conan-release --parallel

r:
	clear
	./build/bin/test_main.exe

clean:
	rm build -fr
