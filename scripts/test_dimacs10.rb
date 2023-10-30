#!/usr/bin/env ruby
DAT_PREFIX = "../data/DIMACS10/"
OPTIONS = "-v --weighted --prob maxcut --max-itr 3000"
P = [0.1, 0.2, 0.3, 0]
SIZE_CAP = 500000
K_CAP = 200

def get_csr_size(fname)
  File.open(fname).each_with_index do |line, i|
    if (i == 2)
      return line.split(/\s+/)[1].to_i
    end
  end
end

Dir.glob(DAT_PREFIX + "*") do |fname|
  next if fname == "." || fname == ".."
  csr_filename = fname
  dat_name = File.basename(fname, "")
  n = get_csr_size(csr_filename)
  kk = [(Math::sqrt(2 * n)).ceil, K_CAP].min

  if n > SIZE_CAP
    print "Test data too large (n = #{n}), skipping... \n"
  else
    P.each do |p|
      if p == 0
        suffix = "_full"
        full = "--full"
      else
        suffix = "_#{p}k"
        full = ""
      end
      real_p = p == 0 ? "" : "#{(kk * p).ceil}"
      print "Testing #{dat_name}... \n"
      cmd = "./rbr #{full} #{OPTIONS} #{csr_filename} #{kk} #{real_p} > logs/#{dat_name}#{suffix}.log"
      %x( #{cmd} )
    end
  end
end
