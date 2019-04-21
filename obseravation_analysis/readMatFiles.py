from meza import io

records = io.read("C:\1UNRuniversityFolder\Dissertation\Chapter 2-snow-forest\obseravation_analysis\Jemez\valles_dep_0810.mat") # only file path, no file objects
print(next(records))