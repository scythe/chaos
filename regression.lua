
require "matrix"

--n-dimensional linear regression
function regression(xt, y)
  --make sure our arrays can be manipulated like matrices
  for k, v in ipairs{xt,y} do _= getmetatable(v) == matrix or v = setmetatable(v, matrix) end
  xtt = transpose(xt)
  beta = inverse(xtt * xt) * xtt * y
  xbar, ybar = {}, 0
  for i = 1, #xt do
    ybar = ybar + y[i]
    for j = 1, #xt[i] do
      xbar[j] = xbar[j] and xbar[j] + xt[i][j] or 0
    end
  end
  y = y / #xt
  for i = 1, #xbar do xbar[i] = xbar[i] / #xt end
  xbar = setmetatable({xbar}, matrix)
  alpha = y - beta * xbar
  return alpha, beta
end

function findplane(xt)
  ind_var_t, dep_var_t = setmetatable({}, matrix), setmetatable({}, matrix)
  for i = 1, #xt do
    ind_var_t[i], dep_var_t[i] = xt[i], {}
    dep_var_t[i] = table.remove(ind_var_t[i]) -- moves the last dimension into the dependent variable table
  end
  return regression(ind_var_t, dep_var_t)
end

function distance(xt, alpha, beta)
  return (xt[1] * beta[1] + xt[2] * beta[2] + xt[3] * beta[3] - xt[4] + alpha) / (beta[1]^2 + beta[2]^2 + beta[3]^2 + 1)^(1/2)
end
