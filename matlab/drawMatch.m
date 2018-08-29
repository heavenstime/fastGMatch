function outImg = drawMatch(inImg, cx, cy, r, thick, theta)

rIO    = floor(r); 
outImg = inImg;
cosTheta = cos(theta);
sinTheta = sin(theta);

for tk = -thick:thick
  rI = rIO + tk;
  for xp = -rI:rI
    yp = floor(sqrt(rI * rI - xp * xp) + 0.01);
    outImg( yp + cy, xp + cx) = 255;
    outImg(-yp + cy, xp + cx) = 255;
  end

  for yp = -rI:rI
    xp = floor(sqrt(rI * rI - yp * yp) + 0.01);
    outImg( yp + cy,  xp + cx) = 255;
    outImg( yp + cy, -xp + cx) = 255;
  end

  for dp = 0:rI
    outImg( floor(dp * sinTheta - cosTheta * tk) + cy, floor(dp * cosTheta + sinTheta * tk) + cx) = 255;
  end
end
  
end
