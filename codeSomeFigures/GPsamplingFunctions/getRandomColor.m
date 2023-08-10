function C=getRandomColor()
  a=1; b=1; c=1;
  while a == 1 && b == 1 && c == 1
     a = round(rand);
     b = round(rand);
     c = round(rand);
  end
  C = [a b c];
end


