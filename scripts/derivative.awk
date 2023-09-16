BEGIN{
	LINES = 0;
}{
	n++;
	if (n == 1) {
		y1 = $2;
		x1 = $1;
	} else if (n == 2) {
		y2 = $2;
		x2 = $1;
		h = x2 - x1;
		print x1, (y2-y1)/h;
	} else if (n == LINES) {
		x = x2;
		y3 = $2;
		print x2, 0.5*(y3-y1)/h;
		
		y1 = y2;
		y2 = $2;
		x = $1;
		
		print x, (y2-y1)/h;
	} else {
		x = x2;
		y3 = $2;
		print x2, 0.5*(y3-y1)/h;
		y1 = y2;
		y2 = y3;
		x2 = $1;
	}
}
