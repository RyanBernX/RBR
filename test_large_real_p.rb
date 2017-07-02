EXEC = "src/ex04"
PREFIX = "../Real_Data"
REGEX_COMM = /nC: (\d+)/
REGEX_TIME = /Time: (\d+\.\d+)/
REGEX_Q = /Q: (\d+\.\d+)/
REGEX_CC = /CC: (\d+\.\d+)/
REGEX_S = /S: (\d+\.\d+)/

DATASET = ['amazon', 'youtube', 'dblp', 'livejournal', 
           'email-Enron', 'loc-brightkite_edges', 'loc-gowalla_edges',
           'delaunay_n24', 'europe_osm', 'road_usa']
#THREADS = [12, 8, 4, 2, 1]
THREADS=[12]
NC = [2000, 3000, 4000]
#P = [1, 2, 5, 10, 20]
P = [5]
#NC = [5, 20, 30, 40]
#DATASET = ['amazon', 'dblp']
#THREADS = [4, 8, 12]

P.each do |vp|
NC.each do |nc|
THREADS.each do |th|
  texfile = "output/large_k#{nc}_th#{th}_p#{vp}.tex"
  logfile = "output/large_k#{nc}_th#{th}_p#{vp}.out"
  f = open(texfile, "w")
  fout = open(logfile, "w")
DATASET.each_with_index do |data, i|
  input = "#{PREFIX}/#{data}"
  command = "taskset -c 0-#{th-1} #{EXEC} #{input} #{nc} #{vp} #{th} --max-itr=75 --output=labels/#{data}_#{nc}_#{vp}.label"
  #command = "mpiexec -np 1 --bind-to core:#{th} --map-by core:#{th} #{EXEC} #{input} #{nc} #{th}"
  ret = `#{command}`
  puts ret
  fout.write ret
  REGEX_COMM =~ ret; nc_r = $1.to_i
  REGEX_TIME =~ ret; time = $1.to_f
  REGEX_Q    =~ ret; q    = $1.to_f
  REGEX_CC   =~ ret; cc   = $1.to_f
  REGEX_S    =~ ret; str  = $1.to_f

  data1 = data.gsub("_", '\_')
  f.write sprintf("%s & %d & %.3f & %.3f & %.3f & %.3f \\\\\n", data1, nc_r, cc, str, q, time)
end
f.close
fout.close
end
end
end
