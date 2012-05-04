
require "ppmlib"
require "matrix"
require "regression"

--[[  performs bogacki-shampine integration of 
      the system of differential equations described by fun_t
      with initial values given by iv_t, from time st to time fi
      returns a table of values describing the trajectory.
      the equations should be functions that take a time 
      and a table of prior values, and return a number.]]--
function bogacki_shampine(eqns, ivs, tol, st, fi, maxvs)
  local vals, kvals = {{}}, {{}, {}, {}, {}}
  local vn = #ivs -- the number of variables. we don't need to save this, but using #ivs is confusing.
  
  for i = 1, vn do 
    vals[1][i] = ivs[i]
  end
  local ts = {st}
  local cnt = 1 -- variable stepsize means it saves time to count steps
  
  function addscaledvectors(y, h, k, ...)
    if not h then return y end
    local res = {}
    for i = 1, #y do res[i] = y[i] + h * k[i] end
    return addscaledvectors(res, ...)
  end
  
  local function sumsq(t)
    res = 0
    for i = 1, #t do res = res + t[i]^2 end
    return res
  end
  
  local h = 0.1
  
  repeat
    local t = ts[cnt]
    local cv = vals[cnt]
    for i = 1, vn do kvals[1][i] = eqns[i](t, cv) end
    for i = 1, vn do kvals[2][i] = eqns[i](t+h/2, addscaledvectors(cv, h/2, kvals[1])) end
    for i = 1, vn do kvals[3][i] = eqns[i](t+3*h/4, addscaledvectors(cv, 3*h/4, kvals[2])) end
    local oest = addscaledvectors(cv, 2*h/9, kvals[1], h/3, kvals[2], 4*h/9, kvals[3])
    for i = 1, vn do kvals[4][i] = eqns[i](t+h, oest) end
    local nest = addscaledvectors(cv, 7*h/24, kvals[1], h/4, kvals[2], h/3, kvals[3], h/8, kvals[4])
    local err = math.abs(math.sqrt(sumsq(addscaledvectors(nest, -1, oest))))

    if err == 0 then
      print("what the fuck")
      print(unpack(nest))
      print(unpack(oest))
      break
    end

    if err < tol then -- it worked!
      --print("successful iteration: " .. h .. " " .. err)
      cnt = cnt + 1
      ts[cnt] = t+h
      vals[cnt] = {}
        -- if a variable is periodic, we let it "wrap around"
      for i = 1, vn do vals[cnt][i] = maxvs[i] and nest[i] % maxvs[i] or nest[i] end 
    else
      print("failed: " .. h .. " " .. err)
    end
    h = 0.8 * h * (tol / err)^(1/4)
  until ts[cnt] > fi
  
  return ts, vals
end
    
function complexlorenz(rho1, rho2, b, sigma, e)
  return {
function(t, xt) return -sigma * (xt[1] - xt[2] * math.cos(xt[3])) end,
function(t, xt) return -xt[2] + (rho1 - xt[4]) * xt[1] * math.cos(xt[3]) end,
function(t, xt) return -e - (sigma*xt[2]/xt[1] + (rho1 - xt[4])*xt[1]/xt[2]) * math.sin(xt[3]) end,
function(t, xt) return -b * xt[4] + xt[1] * xt[2] * math.cos(xt[3]) end
}
end

test_times, test_vals = bogacki_shampine(complexlorenz(28, 0, 8/3, 10, 1/10), {1, 1, 1, 1}, 0.0001, 0, 1250, {nil, nil, math.pi * 2, nil})

cartesian = {}

--project the polar coordinate flow onto cartesian coordinates so it can be graphed
for ind, point in ipairs(test_vals) do
  cartesian[ind] = {}
  cartesian[ind][1] = point[1] * math.cos(point[3])
  cartesian[ind][2] = point[1] * math.sin(point[3])
  cartesian[ind][3] = point[2]
  cartesian[ind][4] = point[4]
end

--now we need to find slices. we choose some points, and divide them at random into two slices

slices = {{}, {}}

for i = 1, 100 do 
  slices[i > 50 and 2 or 1][i] = cartesian[math.random(#cartesian)]
end

alphas, betas = {}, {}

repeat
  switchless = true
  for k, slice in ipairs(slices) do alphas[k], betas[k] = findplane(slice) end
  for k, slice in ipairs(slices) do
    l = 3 - k --the index of the other slice
    for i, point in ipairs(slice) do
      --if a point is closer to the other slice than the slice that it is in, move it to the other slice
      if distance(point, alpha[k], beta[k]) > distance(point, alpha[l], beta[l]) then
        table.insert(slices[l], table.remove(slices[k], i))
              --once there are no points to move, we have the best slices we can find
        switchless = false
      end
    end
  end
until switchless

--now alphas and betas represent the desired slices for the atlas

raw_projections = {{}, {}}

for k, projection in ipairs(raw_projections) do
  for i, point in ipairs(cartesian) do
    for dir, mag in ipairs(point) do
      projection[i][dir] = x - (betas[k] and betas[k][dir] or 1) * distance(x, alphas[k], betas[k]) / sqrt(betas[k][1]^2 + betas[k][2]^2 + betas[k][3]^2 + 1)
    end
  end
end


      

