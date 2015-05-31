dyn.load("cauchyUtils.dll")
is.loaded("testExact")

size = 10;
#location = 0.0;
#scale = 1.0;
iterations_num = 1000000;

Out = .C("testExact", location = as.double(0.0), scale = as.double(1.0), size = as.integer(size), iterations_num = as.integer(iterations_num)); 
Out = .C("testExact", location = as.double(0.0), scale = as.double(5.0), size = as.integer(size), iterations_num = as.integer(iterations_num)); 
Out = .C("testExact", location = as.double(3.0), scale = as.double(1.0), size = as.integer(size), iterations_num = as.integer(iterations_num)); 
Out = .C("testExact", location = as.double(3.0), scale = as.double(5.0), size = as.integer(size), iterations_num = as.integer(iterations_num)); 
#Out

dyn.unload("cauchyUtils.dll")