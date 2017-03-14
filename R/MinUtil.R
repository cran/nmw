Grad = function(func, x)
{
# Numerical Gradient by Kyun-Seop Bae (ksbae@amc.seoul.kr)
# Ref: Richardson Extrapolation
# This is an optimized version for languages like C/C#/C++

  n     = length(x)
  x1    = vector(length=n)
  x2    = vector(length=n)
  ga    = vector(length=4)
  gr    = vector(length=n)

  for (i in 2:n) x1[i] = x2[i] = x[i]

  for (i in 1:n) {
    axi = abs(x[i])
    if (axi <= 1) hi = 1e-4
    else          hi = 1e-4 * axi

    for (k in 1:4) {
      x1[i] = x[i] - hi
      x2[i] = x[i] + hi
      ga[k] = (func(x2) - func(x1)) / (2*hi)
      hi = hi / 2
    }

    ga[1] = (ga[2]*4  - ga[1]) / 3
    ga[2] = (ga[3]*4  - ga[2]) / 3
    ga[3] = (ga[4]*4  - ga[3]) / 3
    ga[1] = (ga[2]*16 - ga[1]) / 15
    ga[2] = (ga[3]*16 - ga[2]) / 15
    gr[i] = (ga[2]*64 - ga[1]) / 63
    x1[i] = x2[i] = x[i]
  }

  return(gr)
}


Hessian = function(func, x)
{
  n  = length(x)
  h0 = vector(length=n)
  x1 = vector(length=n)
  x2 = vector(length=n)
  ha = vector(length=4) # Hessian Approximation
  H  = matrix(NA, nrow=n, ncol=n) # Hessian Matrix

  f0 = func(x)

  for (i in 1:n) {
    x1[i] = x2[i] = x[i]
    axi   = abs(x[i])
    if (axi < 1) h0[i] = 1e-4
    else         h0[i] = 1e-4 * axi
  }

  for (i in 1:n) {
    for (j in i:1) {
      hi = h0[i]
      if (i==j) {
        for (k in 1:4) {
          x1[i] = x[i] - hi
          x2[i] = x[i] + hi
          ha[k] = (func(x1) - 2*f0 + func(x2)) / (hi*hi)
          hi = hi / 2
        }
      } else {
        hj = h0[j]
        for (k in 1:4) {
          x1[i] = x[i] - hi
          x1[j] = x[j] - hj
          x2[i] = x[i] + hi
          x2[j] = x[j] + hj
          ha[k] = (func(x1) - 2*f0 + func(x2) - H[i,i]*hi*hi - H[j,j]*hj*hj) / (2*hi*hj)
          hi = hi / 2
          hj = hj / 2
        }
      }
      w = 4
      for (m in 1:2) {
        for (k in 1:(4-m)) ha[k] = (ha[k+1]*w - ha[k]) / (w - 1)
        w = w * 4
      }
      H[i,j] = (ha[2]*64 - ha[1]) / 63
      if (i != j) H[j,i] = H[i,j] 
      x1[j] = x2[j] = x[j]
    }
    x1[i] = x2[i] = x[i]
  }

  return(H)
}


SqrtInvCov = function(M)
{
  EigenResult = eigen(as.matrix(M))
  EigenVector = EigenResult$vectors
  EigenValues = abs(EigenResult$values)
  return(EigenVector %*% diag(1/sqrt(EigenValues)) %*% t(EigenVector))
}

Mx = function(M, x)
{
  e = eigen(as.matrix(M))
  return(e$vectors %*% diag(abs(e$values)^x) %*% t(e$vectors))
}

# "%^%" = function(S, pow) with(eigen(S), vectors %*% (abs(values)^pow * t(vectors))) 
# S%^%(-0.5)

mat2ltv = function(mat)
{
  return(mat[upper.tri(mat,diag=TRUE)])
}


ltv2mat = function(vec)
{
  LENGTH = length(vec)
  DIM = round((sqrt(8*LENGTH+1)-1)/2,0)
  if (DIM*(DIM+1)/2 != LENGTH) return(NULL)
  mat = matrix(nrow=DIM, ncol=DIM)
  for (m in 1:DIM) {
    for (n in 1:DIM) {
      k = max(m,n)
      l = min(m,n)
      p = k * (k-1) / 2 + l
      mat[m,n] = vec[p]
    }
  }
  return(mat)
}


mat2utv = function(mat)
{
  return(mat[lower.tri(mat,diag=TRUE)])
}


utv2mat = function(vec)
{
  LENGTH = length(vec)
  DIM = as.integer(round((sqrt(8*LENGTH+1)-1)/2,0))
  if (DIM*(DIM+1)/2 != LENGTH) return(NULL)
  mat = matrix(nrow=DIM, ncol=DIM)
  for (m in 1:DIM) {
    for (n in 1:DIM) {
      k = min(m,n)
      l = max(m,n)
      p = (2*DIM - k + 2)*(k - 1)/2 + l - k + 1
      mat[m,n] = vec[p]
    }
  }
  return(mat)
}


ScaleVar = function(VarMat, dim1)
{
  M1 = chol(VarMat)
  V1 = diag(M1)
  M2 = abs(10 * (M1 - diag(V1, nrow=dim1))) + diag(V1/exp(0.1), nrow=dim1)
  return(t(M2))
}


DesclVar = function(mUCP, mSCL)
{
  nRow = dim(mUCP)[1]
  maT = matrix(nrow=nRow, ncol=nRow)

  for (i in 1:nRow) {
    for (j in 1:nRow) {
      if (i==j) {
        maT[i,j] = exp(mUCP[i,j]) * mSCL[i,j]
      } else if(i > j) {
        maT[i,j] = mUCP[i,j] * mSCL[i,j]
      } else {
        maT[i,j] = 0
      }
    }
  }
  return(maT %*% t(maT))
}


DECN = function(UCP)
{
#Extern: nTheta, nEta, nEps
#        ThetaL[nTheta]==LB, ThetaI[nTheta]==IE, ThetaU[nTheta]==UB
#        alpha, OMscl, SGscl
  uTheta = UCP[1:e$nTheta]
  Theta = exp(uTheta - e$alpha)/(exp(uTheta - e$alpha) + 1)*(e$UB - e$LB)
  uvOM  = UCP[(e$nTheta + 1):(e$nTheta + e$nEta*(e$nEta + 1)/2)]
  umOM  = utv2mat(uvOM)
  OM    = DesclVar(umOM, e$OMscl)
  ltvOM = mat2ltv(OM)

  udSG  = UCP[(e$nTheta + e$nEta*(e$nEta + 1)/2 + 1):(e$nTheta + e$nEta*(e$nEta + 1)/2 + e$nEps)]
  umSG  = diag(udSG, nrow=e$nEps)
  SG    = DesclVar(umSG, e$SGscl)
  dgSG  = diag(SG)
  return(c(Theta, ltvOM, dgSG))
}

