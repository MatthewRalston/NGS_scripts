#!/home/mrals/.rvm/rubies/ruby-2.0.0-p247/bin/ruby
def main
  while line=STDIN.gets
    if line[0..2] == "HWI"
      line=line.split("\t")
      STDOUT.puts ([line[0].split("_")[0]]+line[1..line.size]).join("\t")
    else
      STDOUT.puts line
    end
  end
end



main()
