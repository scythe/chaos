
require "ppmlib"

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

test_times, test_vals = bogacki_shampine(complexlorenz(28, 0, 8/3, 10, 1/10), {1, 1, 1, 1}, 0.00001, 0, 6250)

xhat = {0.8, -0.6}
yhat = {0.6, 0.8}
zhat = {0, 1}
image = plot3dtbl(test_vals, 800, 800, 330, 430, 150, 300, {xhat, yhat, zhat})

print(#test_vals)

filename = "lorenz"

imfile = io.open(filename .. ".ppm", "w")
ppmhead(imfile, 800, 800)
ppmwritegraph(imfile, image, {0.2, 0.6, 0.8}, 255)
io.close(imfile)