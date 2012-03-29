
function ppmhead(file, width, height)
  file:write("P6\n" .. width .. "\n" .. height .. "\n255\n")
end

ppmhead(io.stdout, 100, 100)

for i = 1, 100 do
  for j = 1, 100 do
    io.write("\0")
    io.write(string.char(254))
    io.write("\0")
  end
end