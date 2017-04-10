/**********************************
* BoxDivN2��BoxDivN3�����������άMeanFilterʱ��
* ��Ӧ��ÿһ�����꣬����͵ĸ���������������õ�ƽ��ֵʱ�õ�����Χ���ظ���
* outΪ�����width*heightΪ��MeanFilter�Ķ�ά���ݣ�rΪ��ֵ�뾶���������Ĵ��ڴ�СΪ��2*r+1��

* ��Ч���Ͽ���BoxDivN3��N2��΢��һ�㡣���Ƕ���BoxDivN1��2-3����
* N2 N3�����Կ���ʹ��ѭ��չ����SIMD���Ż���
***********************************/

void BoxDivN2(int *out, int width, int height, int r)
{
	int m, n, max;
	max = (2 * r + 1) * (2 * r + 1);   //�����box��������ͼ��Χ�ڵ������
	for (int y = 0, index = -1; y < height; ++y) {
		for (int x = 0; x < width; ++x) {
			if (x >= r && x <= width - r - 1 && y >= r && y <= height - r - 1) {
				out[++index] = max;
			} else {
				m = MINT(r,x, width - 1 - x);
				n = MINT(r,y, height - 1 - y);
				out[++index] = (r + 1 + m) * (r + 1 + n);//���ڱ�Ե��
			}

		}
	}   
}

void BoxDivN3(int *out, int width, int height, int r)
{
	int yCnt, xCnt, index = -1, max = (2*r+1)*(2*r+1);
	for(int y = 0; y < r; ++y)
	{
		yCnt = r + 1 + y;
		for (int x = 0; ; ++x)
		{
			if(x < r)
			{
				xCnt = r + 1 +x;
			}
			else if(x < width -r -1)
			{
				xCnt = 2*r+1;
			}
			else if(x < width)
			{
				xCnt = r + 1 + (width - x - 1);
			}
			else
			{
				break;
			}
			out[++index] = yCnt * xCnt;
		}
	}
	//
	yCnt = 2*r + 1;
	for(int y = r; y < height - r - 1; ++y)
	{
		for (int x = 0; ; ++x)
		{
			if(x < r)
			{
				xCnt = r + 1 +x;
			}
			else if(x < width -r -1)
			{
				xCnt = 2*r+1;
				out[++index] = max;
				continue;
			}
			else if(x < width)
			{
				xCnt = r + 1 + (width - x - 1);
			}
			else
			{
				break;
			}
			out[++index] = yCnt * xCnt;
		}
	}
	for(int y = height - r - 1; y < height; ++y)
	{
		yCnt =  r + 1 + (height - y - 1);
		for (int x = 0; ; ++x)
		{
			if(x < r)
			{
				xCnt = r + 1 +x;
			}
			else if(x < width -r -1)
			{
				xCnt = 2*r+1;
			}
			else if(x < width)
			{
				xCnt = r + 1 + (width - x - 1);
			}
			else
			{
				break;
			}
			out[++index] = yCnt * xCnt;
		}
	}
}


void BoxDivN1(int *out, int width, int height, int r) {
	int * input = (int *)malloc(sizeof(int)*width * height);
	for (int i = 0; i < width * height; ++i)
	{
		input[i] = 1;
	}
	BoxFilter(input, out,r,width,height, 1);
    delete [] input;
}


template<class T1, class T2>
void BoxFilter(T1 *data, T2 *outdata, int r, int width, int height, T2 t) {
	T2 *cum = (T2 *) malloc(sizeof(T2) * width * height);
	if (!cum) {
		//LOG("BoxFilter malloc failed!!!!!!");
		return;
	}
	/**********this  ringht this*************/

	const int len = width * height;
	const int block = 4; //

	//cum y
	for (int i = 0; i < width; ++i) {//the first row
		cum[i] = data[i];
	}
	for (int i = width; i < len; i += width) {
		int end = i + width / block * block;
		for (int j = i; j < end; j += block) {
			cum[j] = cum[j - width] + data[j];
			cum[j + 1] = cum[j - width + 1] + data[j + 1];
			cum[j + 2] = cum[j - width + 2] + data[j + 2];
			cum[j + 3] = cum[j - width + 3] + data[j + 3];
		}
		for (int j = end; j < i + width; ++j) {
			cum[j] = cum[j - width] + data[j];
		}
	}
	//diff y
	const int R_WIDTH = r * width;
	const int R1_WIDTH = width * (r + 1);
	for (int i = 0 ; i < (r + 1) * width; i += block)   //����Ҫ����end�ǲ���block�������������㲻�ǣ������Ĳ���Ҳ���ں������¼�����ȷ��ֵ
	{
		outdata[i] = cum[R_WIDTH+i];
		outdata[i + 1] = cum[R_WIDTH+i + 1];
		outdata[i + 2] = cum[R_WIDTH+i + 2];
		outdata[i + 3] = cum[R_WIDTH+i + 3];
	}
	for (int i = (r + 1) * width; i < (height - r - 1) * width; i += block) {
		outdata[i] = cum[i + R_WIDTH] - cum[i - R1_WIDTH];
		outdata[i + 1] = cum[i + R_WIDTH + 1] - cum[i - R1_WIDTH + 1];
		outdata[i + 2] = cum[i + R_WIDTH + 2] - cum[i - R1_WIDTH + 2];
		outdata[i + 3] = cum[i + R_WIDTH + 3] - cum[i - R1_WIDTH + 3];
	}
	for (int i = height - r - 1; i < height; ++i) {
		int end = width / block * block;
		int outIndex = i * width;
		int topIndex = outIndex - R1_WIDTH;
		int bottomIndex = (height - 1) * width;
		for (int y = 0; y < end; y += block) {
			outdata[outIndex] = cum[bottomIndex] - cum[topIndex];
			outdata[outIndex + 1] = cum[bottomIndex + 1] - cum[topIndex + 1];
			outdata[outIndex + 2] = cum[bottomIndex + 2] - cum[topIndex + 2];
			outdata[outIndex + 3] = cum[bottomIndex + 3] - cum[topIndex + 3];
			outIndex += block;
			topIndex += block;
			bottomIndex += block;
		}
		for (int y = end; y < width; ++y) {
			outdata[outIndex++] = cum[bottomIndex++] - cum[topIndex++];
		}
	}


	//cum x
	for (int y = 0; y < width * height; y += width) {
		cum[y] = outdata[y];  //�����һ��
	}
	for (int y = 0; y < height / 4 * 4; y += 4) {
		//y01234����ÿ�е�����
		int y0 = y * width, y1 = (y + 1) * width, y2 = (y + 2) * width, y3 = (y + 3) * width, y4 =
			(y + 4) * width;
		//ѭ��չ����ÿ�δ������С���ÿ�δӵڶ���Ԫ�ؿ�ʼ��+1������Ϊ��һ���Ѿ��������
		for (int i = y0 + 1; i < y1; ++i) {  //����һ��
			cum[i] = outdata[i] + cum[i - 1];
		}
		for (int i = y1 + 1; i < y2; ++i) {
			cum[i] = outdata[i] + cum[i - 1];
		}
		for (int i = y2 + 1; i < y3; ++i) {
			cum[i] = outdata[i] + cum[i - 1];
		}
		for (int i = y3 + 1; i < y4; ++i) {
			cum[i] = outdata[i] + cum[i - 1];
		}
	}
	for (int y = height / 4 * 4; y < height; ++y) {  //����ѭ��չ�����ʣ����
		for (int i = y * width + 1; i < (y+1)*width; ++i) {
			cum[i] = outdata[i] + cum[i - 1];
		}
	}
	//diff x
	for (int y = 0; y < height; ++y) {
		for (int x = 0; x < r + 1; ++x) {
			outdata[y * width + x] = cum[y * width + x + r];
		}
		for (int x = r + 1; x < width - r; ++x) {
			outdata[y * width + x] =
				(cum[y * width + x + r] - cum[y * width + x - r - 1]);
		}
		for (int x = width - r; x < width; ++x) {
			outdata[y * width + x] = (cum[y * width + width - 1] - cum[y * width + x - r - 1]) ;
		}
	}

	delete [] cum;
	cum = NULL;
}