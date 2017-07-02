EXEC = "taskset -c 0-0 src/ex02"
PREFIX = "../Real_Data"
REGEX_TIME = /Time: (\d+\.\d+)/
REGEX_MIS = /Mis%: (\d+\.\d+)/
REGEX_OBJ = /Obj: (\d+\.\d+[Ee][+-]\d+)/

NC = [8, 2, 4]
DATASET = ['Caltech36', 'polblogs', 'Simmons81']

DATASET.each_with_index do |data, i|
  input = "#{PREFIX}/#{data}"
  info = "#{PREFIX}/#{data}.info.converted"
  command = "#{EXEC} #{input} #{NC[i]} #{info} #{data}_mis labels/#{data}.label"
  ret = `#{command}`
  puts ret
  REGEX_TIME =~ ret; time = $1.to_f
  REGEX_MIS  =~ ret; mis  = $1.to_f
  REGEX_OBJ  =~ ret; obj  = $1.to_f
  #puts sprintf("%s & %.3f & %2.2f & %3.2e \\\\", data, time, mis, obj)
end
