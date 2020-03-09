#include "stdafx.h"
#include <iostream>
#include <opencv2/opencv.hpp>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include "opencv2/core/core.hpp"
#include "opencv2/features2d/features2d.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/calib3d/calib3d.hpp"
using namespace cv;
using namespace std;
//String windowName = "Okarwa"; //Name of the window

//namedWindow(windowName); // Create a window

//imshow(windowName, image); // Show our image inside the created window.

//waitKey(0); // Wait for any keystroke in the window

//destroyWindow(windowName); //destroy the created window


struct Figury {
	Mat figura;
	Figury *next;
};

void tworzymy(Figury *&start, Mat pierwszy, Mat drugi, Mat trzeci, Mat czwarty, Mat piaty, Mat szosty, Mat siodmy, Mat osmy, Mat dziewiaty)
{
	Figury *obecny;
	start = NULL;
	obecny = new Figury;
	obecny->figura =pierwszy;
	obecny->next = start;
	start = obecny;
	obecny = new Figury;
	obecny->figura = drugi;
	obecny->next = start;
	start = obecny;
	obecny = new Figury;
	obecny->figura = trzeci;
	obecny->next = start;
	start = obecny;
	obecny = new Figury;
	obecny->figura = czwarty;
	obecny->next = start;
	start = obecny;
	obecny = new Figury;
	obecny->figura = piaty;
	obecny->next = start;
	start = obecny;
	obecny = new Figury;
	obecny->figura =szosty;
	obecny->next = start;
	start = obecny;
	obecny = new Figury;
	obecny->figura = siodmy;
	obecny->next = start;
	start = obecny;
	obecny = new Figury;
	obecny->figura =osmy;
	obecny->next = start;
	start = obecny;
	obecny = new Figury;
	obecny->figura = dziewiaty;
	obecny->next = start;
	start = obecny;
}


int main(int argc, char** argv)
{
	int erozja = 3;
	Mat elementt = getStructuringElement(MORPH_CROSS,
		Size(2 * erozja + 1, 2 * erozja + 1),
		Point(erozja, erozja));
	Mat czerwo = imread("C:/Users/Hieronim Kaszanka/Desktop/Wyk³ady i ró¿ne/VI SEMESTR/WMA/Projekt/szukane/figury/czerwo.jpg");
	Mat trefl = imread("C:/Users/Hieronim Kaszanka/Desktop/Wyk³ady i ró¿ne/VI SEMESTR/WMA/Projekt/szukane/figury/trefl.jpg");
	Mat dzwonek = imread("C:/Users/Hieronim Kaszanka/Desktop/Wyk³ady i ró¿ne/VI SEMESTR/WMA/Projekt/szukane/figury/dzwonek.jpg");
	Mat wino = imread("C:/Users/Hieronim Kaszanka/Desktop/Wyk³ady i ró¿ne/VI SEMESTR/WMA/Projekt/szukane/figury/wino.jpg");
	Mat krol = imread("C:/Users/Hieronim Kaszanka/Desktop/Wyk³ady i ró¿ne/VI SEMESTR/WMA/Projekt/szukane/figury/krool.jpg");
	Mat jopek = imread("C:/Users/Hieronim Kaszanka/Desktop/Wyk³ady i ró¿ne/VI SEMESTR/WMA/Projekt/szukane/figury/jopek.jpg");
	Mat dama = imread("C:/Users/Hieronim Kaszanka/Desktop/Wyk³ady i ró¿ne/VI SEMESTR/WMA/Projekt/szukane/figury/dama.jpg");
	Mat dyszka = imread("C:/Users/Hieronim Kaszanka/Desktop/Wyk³ady i ró¿ne/VI SEMESTR/WMA/Projekt/szukane/figury/10.jpg");
	Mat as = imread("C:/Users/Hieronim Kaszanka/Desktop/Wyk³ady i ró¿ne/VI SEMESTR/WMA/Projekt/szukane/figury/ass.jpg");
	Figury *start;
	tworzymy(start, as, dyszka, dama, jopek, krol, wino, dzwonek, trefl, czerwo);
	Mat image = imread("C:/Users/Hieronim Kaszanka/Desktop/Wyk³ady i ró¿ne/VI SEMESTR/WMA/Projekt/Zdjecia/3/1.jpg");
	if (image.empty()) // Check for failure
	{
		cout << "Could not open or find the image" << endl;
		system("pause"); //wait for any key press
		return -1;
	}
	else
	{
		cout << "dziala" << endl;
	}
	Size size(1500, 1000);
	resize(image, image, size);//resize image
	String windowNameee = "elo"; //Name of the window
	namedWindow(windowNameee); // Create a window
	imshow(windowNameee, image); // Show our image inside the created window.
	waitKey(0); // Wait for any keystroke in the window
	destroyWindow(windowNameee); //destroy the created window
	Mat mask(image.rows, image.cols, CV_8UC1);
	mask.setTo(255);
	rectangle(mask, Point(0,0), Point(1500, 412), Scalar(0), 500, -1);
	Mat imask= image.clone();
	image.setTo(0, mask);
	for (int y = 0; y < image.rows; y++)
	{
		for (int x = 0; x < image.cols; x++)
		{
			Vec3b piksel = image.at<Vec3b>(Point(x, y)); // 3  bity BGR
			if (!(piksel[2] < 252) || (piksel[1] > 120 && piksel[1] < 252) || (piksel[0] >150 && piksel[0] < 252))
			{
				piksel[0] = 255;
				piksel[1] = 255;
				piksel[2] = 255;
				image.at<Vec3b>(Point(x, y)) = piksel;
			}
		}
	}


	RNG rng(12345);
	int rozmiarerozji = 2 ;
	Mat element = getStructuringElement(MORPH_CROSS,
		Size(2 * rozmiarerozji , 2 * rozmiarerozji),
		Point(rozmiarerozji, rozmiarerozji));

	cvtColor(image, image, CV_BGR2GRAY); //szaro
	threshold(image, image, 220, 255, THRESH_BINARY); //binaryzujemy
	erode(image, image, element);
	medianBlur(image, image, 3);
	Canny(image, image, 50, 50 * 2, 5); // krawêdzie 
	vector<vector<Point> > kontury; //kazdy kontur jest wektorem punktów
	findContours(image, kontury, CV_RETR_TREE, CV_CHAIN_APPROX_NONE); //hierarchia 
	Mat rys = Mat::zeros(image.size(), CV_8UC3); // rys wype³niony zerami-czarny
	int licznik = 1;
	int toto=0;
	for (int i = 0; i < kontury.size(); i++)
	{
		double area = contourArea(kontury[i]); //liczenie obszaru konturu
		if (area > 1000 && area < 5000 && licznik %2==0)
		{
			Scalar color = Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
			drawContours(rys, kontury, i, color, 2, 8); //(image, contours, int contourIdx, const Scalar& color, int thickness=1, int lineType=8, InputArray hierarchy=noArray(), int maxLevel=INT_MAX, Point offset=Point() )	
			toto++;


		}
		licznik++;
	
	} //Dzia³a!!!!!! PogChamp
	int finito = 1000; 
	vector<vector<Point> > pomocny;
	vector<vector<Point> > serce;
	vector<vector<Point> > trefl1;
	vector<vector<Point> > sicko;
	vector<vector<Point> > ding;
	vector<vector<Point> > winko;
	vector<vector<Point> > krul;
	vector<vector<Point> > walet;
	vector<vector<Point> > damulka;
	vector<vector<Point> > ten;
	vector<vector<Point> > ace;
	vector<vector<Point> >hullszukane(1);

	RotatedRect ok;


	for (int j = 0; j < finito; j++)
	{
		Size size(100, 100);
		resize(start->figura, start->figura, size);//resize image
		vector<vector<Point> > beka;
		cvtColor(start->figura, start->figura, CV_BGR2GRAY);
		threshold(start->figura, start->figura, 175, 255, THRESH_BINARY);
		erode(start->figura, start->figura, element);
		Canny(start->figura, start->figura, 5, 5 * 2, 5);
		findContours(start->figura, beka, CV_RETR_TREE, CV_CHAIN_APPROX_NONE);
		Mat plis = Mat::zeros(start->figura.size(), CV_8UC3);
		Scalar color = Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
		pomocny.push_back(beka[0]);
		double obszarr = contourArea(pomocny[0]);
		for (int k = 0; k < beka.size(); k++)
		{
			double obszar = contourArea(beka[k]); //liczenie obszaru konturu
			if (obszar > obszarr)
			{
				pomocny.pop_back();
				pomocny.push_back(beka[k]);
			}
		}
		drawContours(plis, pomocny, -3, color, 2, 8);
		convexHull(Mat(pomocny[0]), hullszukane[0], false);
		drawContours(plis, hullszukane, -3, color, 2, 8);
		vector<vector<Point> >hull(kontury.size());
		licznik = 1;
			Moments M = moments(beka[0], 1);
			double nu11 = M.nu11 * 100;
			double nu20 = M.nu20 * 100;
			double nu02 = M.nu02 * 100;
			double nu21 = M.nu21 * 100;
			double nu12 = M.nu12 * 100;
			double nu03 = M.nu03 * 100;
			double nu30 = M.nu30 * 100;
			double MOM1 = nu20 + nu02;
			double MOM2 = ((nu20 - nu02)*(nu20 - nu02)) + (4 * nu11*nu11);
			double MOM3 = ((nu30 - 3 * nu12)*(nu30 - 3 * nu12)) + ((3 * nu21 - nu03)*(3 * nu21 - nu03));
			double MOM4 = ((nu30 - nu12)*(nu30 - nu12)) + ((nu21 - nu03)*(nu21 - nu03));
			double MOM5 = (nu30 - 3 * nu12)*(nu30 + nu12)*((nu30 + nu12)*(nu30 + nu12) - 3 * (nu21 + nu03)*(nu21 + nu03))+ (3 * nu21 - nu03)*(nu21 + nu03)*(3* (nu30 + nu12)*(nu30 + nu12)- (nu21 + nu03)*(nu21 + nu03));
			double MOM6 = (nu20 - nu02)*((nu30 - nu12)*(nu30 - nu12)- (nu21 + nu03)*(nu21 + nu03))+ 4 * nu11*(nu30 + nu12)*(nu21 + nu03);
			double MOM7 = (nu20*nu02) - (nu11*nu11);


			for (int i = 0; i < kontury.size(); i++)
			{
				double area = contourArea(kontury[i]); //liczenie obszaru konturu
				double wynik;
				double wynik2;
				double perimeter;
				if (area > 1000 && area < 5000 && licznik % 2 == 0)
				{
					convexHull(Mat(kontury[i]), hull[i], false);
					wynik = matchShapes(pomocny[0], kontury[i], CV_CONTOURS_MATCH_I1, 10);
					wynik2 = matchShapes(hullszukane[0], hull[i], CV_CONTOURS_MATCH_I1, 1);
					perimeter = arcLength(kontury[i], true);
					vector<vector<Point> >hull(kontury.size());

					Moments MM = moments(kontury[i], 1);
					double nuu11 = MM.nu11 * 100;
					double nuu20 = MM.nu20 * 100;
					double nuu02 = MM.nu02 * 100;
					double nuu21 = MM.nu21 * 100;
					double nuu12 = MM.nu12 * 100;
					double nuu03 = MM.nu03 * 100;
					double nuu30 = MM.nu30 * 100;
					double M1 = nuu20 + nuu02;
					double M2 = ((nuu20 - nuu02)*(nuu20 - nuu02)) + (4 * nu11*nu11);
					double M3 = ((nuu30 - 3 * nuu12)*(nuu30 - 3 * nuu12)) + ((3 * nuu21 - nuu03)*(3 * nuu21 - nuu03));
					double M4 = ((nuu30 - nuu12)*(nuu30 - nuu12)) + ((nuu21 - nuu03)*(nuu21 - nuu03));
					double M5 = (nuu30 - 3 * nuu12)*(nuu30 + nuu12)*((nuu30 + nuu12)*(nuu30 + nuu12) - 3 * (nuu21 + nuu03)*(nuu21 + nuu03)) + (3 * nuu21 - nuu03)*(nuu21 + nuu03)*(3 * (nuu30 + nuu12)*(nuu30 + nuu12) - (nuu21 + nuu03)*(nuu21 + nuu03));
					double M6 = (nuu20 - nuu02)*((nuu30 - nuu12)*(nuu30 - nuu12) - (nuu21 + nu03)*(nuu21 + nuu03)) + 4 * nuu11*(nuu30 + nuu12)*(nuu21 + nuu03);
					double M7 = (nuu20*nuu02) - (nuu11*nuu11);


					if ((j == 0) && (wynik < 0.15) && (wynik2 < 0.15) && ((perimeter > 100) && (perimeter < 240)) && (M1 > 15) && (M1 < 19) && (M2 > 0) && (M2 <15) && (M3 > 1) && (M3 <15)&& (M4 > 2) && (M4 <7))
					{
						if (area > 1000 && area < 5800)
						{
							serce.push_back(kontury[i]);
							sicko.push_back(kontury[i]);
						}
					}

					if ((j == 1) && (wynik < 0.25) && (wynik2 < 0.25)  && ((perimeter > 180) && (perimeter < 700)) && (M1 > 15) && (M1 < 19) && (M2 > 0) && (M2 <18) && (M3 > 2) && (M3 <11) && (M4 > 0.0) && (M4 <2.4)&&(M7 > 75) && (M7 <85))
					{
						if (area > 1000 && area < 5800)
						{
							trefl1.push_back(kontury[i]);
							sicko.push_back(kontury[i]);
						}
					}

					if ((j == 2) && (wynik < 0.1) && (wynik2 < 0.1) && ((perimeter > 50) && (perimeter < 400)) && (M1 > 0) && (M1 < 30) && (M2 > 0) && (M2 <20) && (M3 >0) && (M3 <2) && (M4 > 0.0) && (M4 <5) && (M7 > 60) && (M7 <80))
					{
						if (area > 1000 && area < 5800)
						{
							ding.push_back(kontury[i]);
							sicko.push_back(kontury[i]);
						}
					}

					if ((j == 3) && (wynik < 0.3) && (wynik2 < 0.3) && ((perimeter > 150) && (perimeter < 250)) && (M1 > 16) && (M1 < 18) && (M2 > 0) && (M2 <12) && (M3 >0) && (M3 <20) && (M4 > 0) && (M4 <2) && (M7 > 70) && (M7 <75))
					{
						if (area > 1000 && area < 5800)
						{
							winko.push_back(kontury[i]);
							sicko.push_back(kontury[i]);
						}
					}


					if ((j == 4) && (wynik < 0.6) && (wynik2 < 0.4) && ((perimeter >320) && (perimeter <600)) && (M1 > 20) && (M1 < 40) && (M2 > 1) && (M2 < 350) && (M3 > 0) && (M3 < 11))
					{
						if (area > 1000 && area < 5800)
						{
							krul.push_back(kontury[i]);
							sicko.push_back(kontury[i]);
						}
					}


					if ((j == 5) && (wynik < 4) && (wynik2 < 4) && ((perimeter >200) && (perimeter <400)) && (M1 > 30) && (M1 < 70) && (M2 > 400) && (M2 < 1500) && (M3 > 160) && (M3 < 4320))
					{
						if (area > 1000 && area < 5800)
						{
							walet.push_back(kontury[i]);
							sicko.push_back(kontury[i]);
						}
					}


					if ((j == 6) && (wynik < 0.3) && (wynik2 < 0.3) && ((perimeter >100) && (perimeter <490)) && (M1 > 18) && (M1 < 21) && (M2 > 12) && (M2 < 157) && (M3 > 1) && (M3 < 5) && (M4 > 0) && (M4 < 3) && (M5 > -1) && (M5 < 1) && (M6 >-4) && (M6 < 1.5) && (M7 > 65) && (M7 < 75) )
					{
						if (area > 1000 && area < 5800)
						{
							damulka.push_back(kontury[i]);
							sicko.push_back(kontury[i]);
						}
					}



					if ((j == 7) && (wynik < 0.15) && (wynik2 < 0.15) && ((perimeter >100) && (perimeter <600)))
					{
						if (area > 1000 && area < 5800)
						{
							ten.push_back(kontury[i]);
							sicko.push_back(kontury[i]);
						}
					}

					if ((j == 8) && (wynik < 0.15) && (wynik2 < 0.15) && ((perimeter >100) && (perimeter <600)))
					{
						if (area > 1000 && area < 5800)
						{
							ace.push_back(kontury[i]);
							sicko.push_back(kontury[i]);
						}
					}

				}
			licznik++;
		}

		if (j == 8)
		{

			Mat help = Mat::zeros(image.size(), CV_8UC3);
			drawContours(help,sicko, -3, color, 2, 8);
			String windowNam = "Okawa"; //Name of the window
			namedWindow(windowNam); // Create a window
			imshow(windowNam, help); // Show our image inside the created window.
			waitKey(0); // Wait for any keystroke in the window
			destroyWindow(windowNam); //destroy the created window
		}

		
	
		pomocny.pop_back();
		start = start->next;
		if (start == NULL)
		{
			finito = j + 1;
		}
	}

	double cx;
	double cy;
	double dx;
	double dy;
	int krolczerwo = 0;
	int kroltrefl = 0;
	int krolding = 0;
	int krolwinko = 0;
	int waletczerwo = 0;
	int walettrefl = 0;
	int waletding = 0;
	int damaczerwo = 0;
	int damawinko = 0;
	int dyszkaczerwo = 0;
	int dyszkading = 0;
	int asczerwo = 0;
	int asding = 0;
	int aswinko = 0;


	for (int z=0; z < krul.size(); z++) {
		Moments mu = moments(krul[z]);
		cx = mu.m10 / mu.m00;
		cy = mu.m01 / mu.m00;
		for (int zz = 0; zz < serce.size(); zz++) {
			Moments muu = moments(serce[zz]);
			dx = muu.m10 / muu.m00;
			dy = muu.m01 / muu.m00;
			double odl = sqrt(abs((cx - dx)*(cx - dx)) + abs((cy - dy)*(cy - dy)));
			if (odl < 100 && krolczerwo < 1 && odl != 0) {
				krolczerwo++;
			}
			
		}
		for (int zz = 0; zz < trefl1.size(); zz++) {
			Moments muu = moments(trefl1[zz]);
			dx = muu.m10 / muu.m00;
			dy = muu.m01 / muu.m00;
			double odl = sqrt(abs((cx - dx)*(cx - dx)) + abs((cy - dy)*(cy - dy)));
			if (odl < 100 && kroltrefl < 1 && odl != 0) {
				kroltrefl++;
			}
			
		}

		for (int zz = 0; zz < ding.size(); zz++) {
			Moments muu = moments(ding[zz]);
			dx = muu.m10 / muu.m00;
			dy = muu.m01 / muu.m00;
			double odl = sqrt(abs((cx - dx)*(cx - dx)) + abs((cy - dy)*(cy - dy)));
			if (odl < 100 && krolding < 1 && odl != 0) {
				krolding++;
			}
		}

		for (int zz = 0; zz < winko.size(); zz++) {
			Moments muu = moments(winko[zz]);
			dx = muu.m10 / muu.m00;
			dy = muu.m01 / muu.m00;
			double odl = sqrt(abs((cx - dx)*(cx - dx)) + abs((cy - dy)*(cy - dy)));
			if (odl < 100 && krolwinko < 1 && odl != 0) {
				krolwinko++;
			}
		}
	}


	for (int z = 0; z < walet.size(); z++) {
		Moments mu = moments(walet[z]);
		cx = mu.m10 / mu.m00;
		cy = mu.m01 / mu.m00;
		for (int zz = 0; zz < serce.size(); zz++) {
			Moments muu = moments(serce[zz]);
			dx = muu.m10 / muu.m00;
			dy = muu.m01 / muu.m00;
			double odl = sqrt(abs((cx - dx)*(cx - dx)) + abs((cy - dy)*(cy - dy)));
			if (odl < 100 && waletczerwo < 1 && odl != 0) {
				waletczerwo++;
			}

		}
		for (int zz = 0; zz < trefl1.size(); zz++) {
			Moments muu = moments(trefl1[zz]);
			dx = muu.m10 / muu.m00;
			dy = muu.m01 / muu.m00;
			double odl = sqrt(abs((cx - dx)*(cx - dx)) + abs((cy - dy)*(cy - dy)));
			if (odl < 115 && walettrefl < 1 && odl != 0) {
				walettrefl++;
			}

		}

		for (int zz = 0; zz < ding.size(); zz++) {
			Moments muu = moments(ding[zz]);
			dx = muu.m10 / muu.m00;
			dy = muu.m01 / muu.m00;
			double odl = sqrt(abs((cx - dx)*(cx - dx)) + abs((cy - dy)*(cy - dy)));
			if (odl < 110 && waletding < 1 && odl != 0) {
				waletding++;
			}
		}

	}


		for (int z = 0; z < damulka.size(); z++) {
			Moments mu = moments(damulka[z]);
			cx = mu.m10 / mu.m00;
			cy = mu.m01 / mu.m00;
			for (int zz = 0; zz < serce.size(); zz++) {
				Moments muu = moments(serce[zz]);
				dx = muu.m10 / muu.m00;
				dy = muu.m01 / muu.m00;
				double odl = sqrt(abs((cx - dx)*(cx - dx)) + abs((cy - dy)*(cy - dy)));
				if (odl < 100 && damaczerwo < 1 && odl != 0) {
					damaczerwo++;
				}

			}
			
			for (int zz = 0; zz < winko.size(); zz++) {
				Moments muu = moments(winko[zz]);
				dx = muu.m10 / muu.m00;
				dy = muu.m01 / muu.m00;
				double odl = sqrt(abs((cx - dx)*(cx - dx)) + abs((cy - dy)*(cy - dy)));
				if (odl < 100 && damawinko < 1 && odl != 0) {
					damawinko++;
				}
			}
		}


		for (int z = 0; z < ten.size(); z++) {
			Moments mu = moments(ten[z]);
			cx = mu.m10 / mu.m00;
			cy = mu.m01 / mu.m00;
			for (int zz = 0; zz < serce.size(); zz++) {
				Moments muu = moments(serce[zz]);
				dx = muu.m10 / muu.m00;
				dy = muu.m01 / muu.m00;
				double odl = sqrt(abs((cx - dx)*(cx - dx)) + abs((cy - dy)*(cy - dy)));
				if (odl < 100 && dyszkaczerwo < 1 && odl != 0) {
					dyszkaczerwo++;
				}

			}
			for (int zz = 0; zz < ding.size(); zz++) {
				Moments muu = moments(ding[zz]);
				dx = muu.m10 / muu.m00;
				dy = muu.m01 / muu.m00;
				double odl = sqrt(abs((cx - dx)*(cx - dx)) + abs((cy - dy)*(cy - dy)));
				if (odl < 100 && dyszkading < 1 && odl != 0) {
					dyszkading++;
				}
			}
		}



		for (int z = 0; z < ace.size(); z++) {
			Moments mu = moments(ace[z]);
			cx = mu.m10 / mu.m00;
			cy = mu.m01 / mu.m00;
			for (int zz = 0; zz < serce.size(); zz++) {
				Moments muu = moments(serce[zz]);
				dx = muu.m10 / muu.m00;
				dy = muu.m01 / muu.m00;
				double odl = sqrt(abs((cx - dx)*(cx - dx)) + abs((cy - dy)*(cy - dy)));
				if (odl < 100 && asczerwo < 1 && odl != 0) {
					asczerwo++;
				}

			}
			
			for (int zz = 0; zz < ding.size(); zz++) {
				Moments muu = moments(ding[zz]);
				dx = muu.m10 / muu.m00;
				dy = muu.m01 / muu.m00;
				double odl = sqrt(abs((cx - dx)*(cx - dx)) + abs((cy - dy)*(cy - dy)));
				if (odl < 100 && asding < 1 && odl != 0) {
					asding++;
				}
			}

			for (int zz = 0; zz < winko.size(); zz++) {
				Moments muu = moments(winko[zz]);
				dx = muu.m10 / muu.m00;
				dy = muu.m01 / muu.m00;
				double odl = sqrt(abs((cx - dx)*(cx - dx)) + abs((cy - dy)*(cy - dy)));
				if (odl < 100 && aswinko < 1 && odl != 0) {
					aswinko++;
				}
			}
		}

		cout << "Rozpoznane karty: " << endl;

		if (krolczerwo == 1) {
			cout << "Krol czerwo" << endl;
		}

		if (krolwinko == 1) {
			cout << "Krol wino" << endl;

		}
		if (krolding == 1) {
				cout << "Krol dzwonek" << endl;
			}

		if (kroltrefl == 1) {
				cout << "Krol trefl" << endl;
			}

		if (waletczerwo == 1) {
			cout << "Walet czerwo" << endl;
		}

		if (waletding == 1) {
			cout << "Walet dzwonek" << endl;
		}

		if (walettrefl == 1) {
			cout << "Walet trefl" << endl;
		}

		if (damaczerwo == 1) {
			cout << "Dama czerwo" << endl;
		}

		if (damawinko == 1) {
			cout << "Dama wino" << endl;

		}

		if (dyszkaczerwo == 1) {
			cout << "10 czerwo" << endl;
		}

		if (dyszkading == 1) {
			cout << "10 dzwonek" << endl;
		}

		if (asczerwo == 1) {
			cout << "AS czerwo" << endl;
		}

		if (aswinko == 1) {
			cout << "AS wino" << endl;

		}
		if (asding == 1) {
			cout << "AS dzwonek" << endl;
		}


		cout << "Uklad kart: " << endl;

		if (krolczerwo == 1 && krolding == 1 && kroltrefl == 1 && krolwinko == 1) {
			cout << "Czworeczka" << endl;
		}

		if (krolczerwo == 1 && waletczerwo == 1 && damaczerwo == 1 && asczerwo == 1 && dyszkaczerwo==1) {
			cout << "Krolewski Poker" << endl;
		}

		if (damaczerwo == 1 && waletczerwo == 1 && damawinko == 1 && waletding == 1 && walettrefl==1) {
			cout << "FULL (3 Walety 2 Damy)" << endl;
		}

		if (krolczerwo == 1 && kroltrefl == 1 && aswinko == 1 && asczerwo == 1 && asding==1) {
			cout << "FULL(3 Asy 2 Krole)" << endl;
		}

		if (krolding == 1 && kroltrefl == 1 && damaczerwo == 1 && waletczerwo == 1 && walettrefl == 1) {
			cout << "2 Pary  ( 2 Walety i 2 Krole)" << endl;
		}

	string d;
	cin >> d;
	return 0;


}