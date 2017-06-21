function pdf=estPDF(mask)
h=1/49*ones(7);
pdf=conv2(mask,h,'same');
pdf=pdf+eps;
end