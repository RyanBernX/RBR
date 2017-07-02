EXEC = "src/ex01"
PREFIX = "../"
REGEX_TIME = /Time: (\d+\.\d+)/
REGEX_MIS = /Mis%: (\d+\.\d+)/

NC = [2, 2, 3, 4]
N = [200, 450, 200, 200]

PGROUPS = [[0.05, 0.10, 0.15, 0.20],
           [0.05, 0.10, 0.15, 0.20],
           [0.15, 0.25, 0.35, 0.45],
           [0.10, 0.20, 0.30, 0.40]]

SHAPE = [1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]

ROUNDS = 5

(0...N.size).each do |i|
  param_p = PGROUPS[i]
  rst_avg = Array.new(param_p.size){ Array.new(SHAPE.size){0} }
  rst_time = Array.new(param_p.size){ Array.new(SHAPE.size){0} }
  nc = NC[i]; n = N[i]
  param_p.each_with_index do |p, ip|
    SHAPE.each_with_index do |s, is|
      ROUNDS.times do |j|
        prefix = sprintf("%sSynthetic_n%d_k%d_csr", PREFIX, nc, n)
        filename = "rounds#{j+1}_p#{p}_shape#{s}"
        puts "processing (#{nc}, #{n}), #{filename}"
        command = "#{EXEC} #{prefix}/#{filename} #{nc}"
        ret = `#{command}`
        #puts ret
        REGEX_MIS =~ ret
        rst_avg[ip][is] += $1.to_f
	REGEX_TIME =~ ret 
        rst_time[ip][is] += $1.to_f
      end
      rst_avg[ip][is] /= ROUNDS
      rst_time[ip][is] /= ROUNDS
    end
  end
  out_avg = "#{PREFIX}rst_#{nc}_#{n}_roundC_12_avg.csv"
  out_time = "#{PREFIX}rst_#{nc}_#{n}_roundC_12_time.csv"
  File.open(out_avg, "w") do |f|
    rst_avg.each do |line|
      f.write line.map{|e| e.to_s}.join(',') + "\n"
    end
  end
  File.open(out_time, "w") do |f|
    rst_time.each do |line|
      f.write line.map{|e| e.to_s}.join(',') + "\n"
    end
  end
end

