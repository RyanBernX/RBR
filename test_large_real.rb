EXEC = "src/ex03"
PREFIX = "../Real_Data"
REGEX_TIME = /Time: (\d+\.\d+)/
REGEX_Q = /Q: (\d+\.\d+)/
REGEX_CC = /CC: (\d+\.\d+)/
REGEX_S = /S: (\d+\.\d+)/

#NC = [100, 200, 300, 400]
DATASET = ['amazon', 'youtube', 'dblp', 'livejournal', \
           'email-Enron', 'loc-brightkite_edges', 'loc-gowalla_edges']
#THREADS = [12, 8, 4, 2, 1]
THREADS=[12]
#NC = [100, 300, 400]
NC = [1000]
#DATASET = ['amazon', 'dblp']
#THREADS = [4, 8, 12]

NC.each do |nc|
THREADS.each do |th|
  texfile = "large_k#{nc}_th#{th}.tex"
  logfile = "large_k#{nc}_th#{th}.out"
  f = open(texfile, "w")
  fout = open(logfile, "w")
DATASET.each_with_index do |data, i|
  input = "#{PREFIX}/#{data}"
  command = "taskset -c 0-#{th-1} #{EXEC} #{input} #{nc} #{th}"
  #command = "mpiexec -np 1 --bind-to core:#{th} --map-by core:#{th} #{EXEC} #{input} #{nc} #{th}"
  ret = `#{command}`
  puts ret
  fout.write ret
  REGEX_TIME =~ ret; time = $1.to_f
  REGEX_Q    =~ ret; q    = $1.to_f
  REGEX_CC   =~ ret; cc   = $1.to_f
  REGEX_S    =~ ret; str  = $1.to_f

  data1 = data.gsub("_", '\_')
  f.write sprintf("%s & %.3f & %.3f & %.3f & %.3f \\\\\n", data1, cc, str, q, time)
end
f.close
fout.close
end
end

