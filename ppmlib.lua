
function ppmhead(file, width, height)
  file:write("P6\n" .. width .. "\n" .. height .. "\n255\n")
end

--writes a colored graph. decmagnum is the magic number in the image matrix for axes, etc, for clarity.
function ppmwritegraph(file, image, color, decmagnum)
  local cfl = function(num) return string.char(math.floor(num)) end
  for i, v in ipairs(image) do 
    if v ~= decmagnum then
      file:write(cfl(v*color[1]) .. cfl(v*color[2]) .. cfl(v*color[3])) 
    else
      file:write(cfl(255) .. cfl(255) .. cfl(255))
    end
  end
end

--creates a decorated graph w/axes. getting the right, xoff, yoff requires experimentation at this point.
function plot3dtbl(valt, width, height, xoff, yoff, axisoriginx, axisoriginy, hats)
  
  local image = {} --create a blank image
  for y = 0, height - 1 do for x = 1, width do image[y*width + x] = 0 end end
  
  for k, v in ipairs(valt) do
    local realx, realy = xoff, yoff
  --  print(unpack(v))
    for dir, mag in ipairs(v) do
      realx, realy = realx + hats[dir][1] * mag * 20, realy + hats[dir][2] * mag * 20
    end
    realy = height - realy
  --  print(realx, realy)
    local coord = math.floor(realy) * width + math.floor(realx)
    if image[coord] and k > 100 and image[coord] < 240 then image[coord] = image[coord] + 80 end
  end
  
  for k, v in ipairs(hats) do
    for i = 1, width * 4/3 do
      local coord = math.floor(height - axisoriginy - i * v[2] / 2) * width + math.floor(axisoriginx + i * v[1] / 2)
      if image[coord] then image[coord] = 255 end
    end
  end  
  return image
end
