function bluetan = colormap_bluetan()

darkness = 0.4;

color1 = darkness*[0,    0,  1];
color2 = darkness*[1, 0.,  0.];

white = [1,1,1];

ncolor = 255;

nmiddle = round( (ncolor-1)/2 + 1 );
bluetan = zeros(ncolor,3);

weight2 = linspace(0,1,nmiddle)';
weight1 = 1 - weight2;

for i = 1:nmiddle
  bluetan(i,:) = weight1(i)*color1 + weight2(i)*white;
end
for iFromEnd = 1:nmiddle
  i = ncolor + 1 - iFromEnd;
  bluetan(i,:) = weight1(iFromEnd)*color2 + weight2(iFromEnd)*white;
end
