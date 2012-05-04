--Classical symmetry-reduced trajectories of the dihydrogen cation.
--We introduce alpha (a), a dimensionless parameter that describes the scaling of the z-axis angular momentum
--with respect to the separation of the protons

require "ppmlib"

function dhc_iter(old, step, a)
  local s = old.s + step*old.vs
  local z = old.z + step*old.vz
  local as = 2/(a*s^3) - a*s * ( 1/(a^2*s^2 + (z+1)^2)^(3/2) + 1/(a^2*s^2 + (z-1)^2)^(3/2) )
  local az = -(z+1)/(a^2*s^2 + (z+1)^2)^(3/2) - (z-1)/(a^2*s^2 + (z-1)^2)^(3/2)
  local vs = old.vs + step*old.as
  local vz = old.vz + step*old.az
  local t = old.t + step
  --print(s, z)
  return {s = s, z = z, vs = vs, vz = vz, as = as, az = az, t = t}
end

function dihydrogencation_leapfrog(s_0, z_0, t_max, step, a)
  local val_t = {{s = s_0, z = z_0, vs = 0, vz = 0, as = 0, az = 0, t = 0}}
  local s, z, vs, vz, as, az
  
  for i = 1, t_max / step do
    val_t[#val_t+1] = dhc_iter(val_t[#val_t], step, a)
  end
  
  return val_t
end

--poincare section at the mean value of the azimuthal radius.
--binary symbolic dynamics is just the sign of z.
function dihydrogen_poincare(a)
  local trajectory = dihydrogencation_leapfrog(4, 4, 5000, 0.001, 1)
  
  sm = 0
  for k, v in ipairs(trajectory) do sm = sm + v.s end
  sm = sm / #trajectory
  
  for k, v in ipairs(trajectory) do
    if k < #trajectory and v.s < sm and trajectory[k+1].s > sm then
      print(v.z < 0 and 0 or 1)
    end
  end
end

dihydrogen_poincare(1)

--[[
plot = {}

for k, v in ipairs(trajectory) do plot[k] = {v.s, v.z} end

print(#plot, unpack(plot[1]))

image = plot3dtbl(plot, 800, 800, 200, 400, 200, 400, {{1, 0}, {0, 1}})

imfile = io.open("dihydrogencation.ppm", "w")
ppmhead(imfile, 800, 800)
ppmwritegraph(imfile, image, {0.3, 0.7, 1.0}, 255)
io.close(imfile)

os.execute("pnmtopng dihydrogencation.ppm > dihydrogencation.png")]]