package main

import (
	"fmt"
	"image"
	"image/color"
	"image/jpeg"
	"io"
	"log"
	"math"
	"os"
	"strconv"
	"strings"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/font"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
)

const (
	pi      float64 = 3.1415
	scale   float64 = 1.41421356
	i_scale         = 1 / scale
)

var Path string
var Dk = map[int][]float64{
	2:  {1.0, 1.0},
	4:  {0.6830127, 1.1830127, 0.3169873, -0.1830127},
	6:  {0.47046721, 1.14111692, 0.650365, -0.19093442, -0.12083221, 0.0498175},
	8:  {0.32580343, 1.01094572, 0.89220014, -0.03957503, -0.26450717, 0.0436163, 0.0465036, -0.01498699},
	10: {0.22641898, 0.85394354, 1.02432694, 0.19576696, -0.34265671, -0.04560113, 0.10970265, -0.00882680, -0.01779187, 4.71742793e-3},
	12: {0.15774243, 0.69950381, 1.06226376, 0.44583132, -0.3199866, -0.18351806, 0.13788809, 0.03892321, -0.04466375, 0.00078325115, 0.00675606236, -0.00152353381},
	14: {0.11009943, 0.56079128, 1.03114849, 0.66437248, -0.20351382, -0.31683501, 0.1008467, 0.11400345, -0.05378245, -0.02343994, 0.01774979, 6.07514995e-4, -2.54790472e-3, 5.00226853e-4},
	16: {0.07695562, 0.44246725, 0.95548615, 0.82781653, -0.02238574, -0.40165863, 6.68194092e-4, 0.18207636, -0.02456390, -0.06235021, 0.01977216, 0.01236884, -6.88771926e-3, -5.54004549e-4, 9.55229711e-4, -1.66137261e-4},
	18: {0.05385035, 0.34483430, 0.85534906, 0.92954571, 0.18836955, -0.41475176, -0.13695355, 0.21006834, 0.043452675, -0.09564726, 3.54892813e-4, 0.03162417, -6.67962023e-3, -6.05496058e-3, 2.61296728e-3, 3.25814671e-4, -3.56329759e-4, 5.5645514e-5},
	20: {0.03771716, 0.26612218, 0.74557507, 0.97362811, 0.39763774, -0.35333620, -0.27710988, 0.18012745, 0.13160299, -0.10096657, -0.04165925, 0.04696981, 5.10043697e-3, -0.01517900, 1.97332536e-3, 2.81768659e-3, -9.69947840e-4, -1.64709006e-4, 1.32354367e-4, -1.875841e-5},
}

type Martix struct {
	mat  [][]float64
	n, m int
}
type coord struct {
	x, y float64
}

func New(n, m int) *Martix {
	var M = Martix{}
	M.n = n
	M.m = m
	M.mat = make([][]float64, n)
	for i := 0; i < n; i++ {
		M.mat[i] = make([]float64, m)
	}
	return &M
}

func Det(mat [][]float64) float64 {
	if len(mat) == 1 {
		return mat[0][0]
	}
	sn, d := 1, float64(0)
	for i, val := range mat[0] {
		d += float64(sn) * val * Det(Ex(mat[1:], i))
		sn *= -1
	}
	return d
}

func Ex(mat [][]float64, key int) [][]float64 {
	res := make([][]float64, len(mat))
	size := len(mat[0]) - 1
	for j, r := range mat {
		buf := make([]float64, size)
		copy(buf[:key], r[:key])
		copy(buf[key:], r[key+1:])
		res[j] = buf
	}
	return res
}
func (M *Martix) Mulconst(constValue float64) {
	for i := 0; i < M.n; i++ {
		for j := 0; j < M.m; j++ {
			M.mat[i][j] *= constValue
		}
	}
}

func (M Martix) Mul(second *Martix) *Martix {
	if M.m != second.n {
		fmt.Fprintln(os.Stderr, "Разный размер матриц")
		return nil
	}
	res := Martix{}
	res.n = M.n
	res.m = second.m
	res.mat = make([][]float64, res.n)
	for i := 0; i < res.n; i++ {
		res.mat[i] = make([]float64, second.m)
		for j := 0; j < second.m; j++ {
			res.mat[i][j] = 0
			for k := 0; k < M.m; k++ {
				res.mat[i][j] += M.mat[i][k] * second.mat[k][j]
			}
		}
	}
	return &res
}

func (M Martix) TMat() *Martix {
	res := New(M.m, M.n)
	for i := 0; i < M.n; i++ {
		for j := 0; j < M.m; j++ {
			res.mat[j][i] = M.mat[i][j]
		}
	}
	return res
}

func (M *Martix) Update(n, m int, value float64) {
	M.mat[n][m] = value
}

func (M *Martix) UpdateArray(value []float64) {
	if M.m*M.n >= len(value) {
		count := 0
		for i := 0; i < M.n; i++ {
			for j := 0; j < M.m; j++ {
				M.mat[i][j] = value[count]
				count++
				if count == len(value) {
					return
				}
			}
		}
	}
}

func (M Martix) Print() {
	for _, slice1 := range M.mat {
		for _, value := range slice1 {
			fmt.Print(value, " ")
		}
		fmt.Println()
	}
}

func (M *Martix) Reverse() *Martix {
	res := New(M.n, M.m)
	for i := 0; i < M.n; i++ {
		for j := 0; j < M.m; j++ {
			res.mat[i][j] = M.AlgSum(i, j)
		}
	}
	fmt.Println("here")
	res = res.TMat()
	res.Mulconst(1.0 / Det(M.mat))
	return res
}

func (M *Martix) Minor(n, m int) *Martix {
	res := New(M.n-1, M.m-1)
	fl1, fl2 := 0, 0
	for i := 0; i < M.n; i++ {
		if i == n {
			fl1 = 1
			continue
		}
		for j := 0; j < M.m; j++ {
			if j == m {
				fl2 = 1
				continue
			}
			res.mat[i-fl1][j-fl2] = M.mat[i][j]
		}
		fl2 = 0
	}
	return res
}

func (M *Martix) AlgSum(n, m int) float64 {
	return math.Pow(-1, float64(n+m)) * Det(M.Minor(n, m).mat)
}

func main() {
START:
	fmt.Println("Режим работы:\n1-Сигнал\n2-Изображение\n0-Выход")
	var switcher int
	fmt.Scan(&switcher)
	switch switcher {
	case 1:
		os.MkdirAll("plots/", 0777)
		var ch int
		fmt.Println("1.Статический сигнал\n2.Нестатический")
		fmt.Scan(&ch)
		switch ch {
		case 1:
			os.MkdirAll("plots/Stat/", 0777)
			//Работа с сигналом
			fmt.Print("Размер сигнала - ")
			var size float64
			fmt.Scan(&size)
			fmt.Print("Количество частей сигнала - ")
			var parts float64
			fmt.Scan(&parts)
			name := SignalCreateStat(size, parts)
			signal, err := OpenFile(name, int(parts))
			if err != nil {
				log.Fatal((err))
				os.Exit(1)
			}
			var Yn []float64
			var Xn []float64
			for _, val := range *signal {
				Yn = append(Yn, val.y)
				Xn = append(Xn, val.x)
			}
			fmt.Println("Выберите количество фильтров в вейвлете Добеши(2,4,6,8,10,12,14,16,18,20)")
			var dlvl int
			fmt.Scan(&dlvl)
			for i := 1; ; i++ {
				_, dob := DaubechiNormalize(Yn, dlvl, i)
				var take int = len(dob) / int(math.Pow(2, float64(i)))
				if dob == nil {
					break
				}
				Path = "plots/Stat/Pic_" + strconv.Itoa(i) + ".png"
				plot := CreatePlot("Стационарный сигнал (уровень "+strconv.Itoa(i)+") D"+strconv.Itoa(dlvl), 1000, 600)
				AddToPlot(Xn, Yn, plot, color.RGBA{R: 0, G: 0, B: 0, A: 255})
				AddToPlot(Xn[:take], dob[:take], plot, color.RGBA{R: 0, G: 0, B: 255, A: 255})
				AddToPlot(Xn[take:], dob[take:], plot, color.RGBA{R: 255, G: 0, B: 0, A: 255})

				dob = append(dob[:take], make([]float64, len(dob)-take)...)

				recover := DaubechiRecoverNormalize(dob, dlvl, i)
				Path = "plots/Stat/Pic_Rec" + strconv.Itoa(i) + ".png"
				plot1 := CreatePlot("Восстановленный сигнал (уровень "+strconv.Itoa(i)+") D"+strconv.Itoa(dlvl), 1000, 600)
				AddToPlot(Xn, Yn, plot1, color.RGBA{R: 0, G: 0, B: 0, A: 255})
				AddToPlot(Xn, recover, plot1, color.RGBA{R: 0, G: 255, B: 0, A: 255})
			}
		case 2:
			os.MkdirAll("plots/notStat/", 0777)
			fmt.Print("Размер сигнала - ")
			var size float64
			fmt.Scan(&size)
			fmt.Print("Количество частей сигнала - ")
			var parts float64
			fmt.Scan(&parts)
			name := SignalCreate(size, parts)
			signal, err := OpenFile(name, int(parts))
			if err != nil {
				log.Fatal((err))
				os.Exit(1)
			}
			var Yn []float64
			var Xn []float64
			for _, val := range *signal {
				Yn = append(Yn, val.y)
				Xn = append(Xn, val.x)
			}
			fmt.Println("Выберите количество фильтров в вейвлете Добеши(2,4,6,8,10,12,14,16,18,20)")
			var dlvl int
			fmt.Scan(&dlvl)
			for i := 1; ; i++ {
				_, dob := DaubechiNormalize(Yn, dlvl, i)
				var take int = len(dob) / int(math.Pow(2, float64(i)))
				fmt.Println(i, "=", len(dob))
				if dob == nil {
					break
				}

				Path = "plots/NotStat/Pic_" + strconv.Itoa(i) + ".png"
				plot := CreatePlot("Нестационарный сигнал (уровень "+strconv.Itoa(i)+") D"+strconv.Itoa(dlvl), 1000, 600)
				AddToPlot(Xn, Yn, plot, color.RGBA{R: 0, G: 1, B: 0, A: 255})
				AddToPlot(Xn[:take], dob[:take], plot, color.RGBA{R: 0, G: 0, B: 255, A: 255})
				AddToPlot(Xn[take:], dob[take:], plot, color.RGBA{R: 255, G: 0, B: 0, A: 255})

				dob = append(dob[:take], make([]float64, len(dob)-take)...)

				recover := DaubechiRecoverNormalize(dob, dlvl, i)
				Path = "plots/NotStat/Pic_Rec" + strconv.Itoa(i) + ".png"
				plot1 := CreatePlot("Восстановленный сигнал (уровень "+strconv.Itoa(i)+") D"+strconv.Itoa(dlvl), 1000, 600)
				AddToPlot(Xn, Yn, plot1, color.RGBA{R: 0, G: 0, B: 0, A: 255})
				AddToPlot(Xn, recover, plot1, color.RGBA{R: 0, G: 255, B: 0, A: 255})
			}
		}

	case 2:
		//Работа с картикной здесь
		os.MkdirAll("pics/", 0777)
		picname := "pic8"
		picadress := "pics/" + picname + ".jpg"

		dataPic, size := GetDataPic_1(picadress)
		DataMatrixRed := New(size.X, size.Y)
		for i := 0; i < size.Y; i++ {
			DataMatrixRed.UpdateArray(dataPic[0])
		}
		MatrixPicDataR := New(size.X, size.Y)
		fmt.Println("Выберите вейлвет\n1:Хаара\n2:Добеши\n3:Хаара без растяжения\n4:Добеши без растяжения\n5:Восстановление")
		var flag int
		method := "Haara"
		fmt.Scan(&flag)
		var Db int
		Db = 2
		if flag == 2 {
			fmt.Println("Выберите количество фильтров в вейвлете Добеши(2,4,6,8,10,12,14,16,18,20)")
			fmt.Scan(&Db)
			method = "Daubechi_" + strconv.Itoa(Db)
		}
		if flag == 4 {
			fmt.Println("Выберите количество фильтров в вейвлете Добеши(2,4,6,8,10,12,14,16,18,20)")
			fmt.Scan(&Db)
			os.MkdirAll("pics/"+picname+"noStr", 0777)
			method = "Daubechi_" + strconv.Itoa(Db)
			os.MkdirAll("pics/"+picname+"noStr"+"/"+method, 0777)
		}
		if flag == 5 {
			fmt.Println("Выберите количество фильтров в вейвлете Добеши(2,4,6,8,10,12,14,16,18,20)")
			fmt.Scan(&Db)
			fmt.Println("Уровень сжатия картинки")
			var startlvl int
			fmt.Scan(&startlvl)
			imgfile, err := os.Open("pics/" + picname + "noStr/Daubechi_" + strconv.Itoa(Db) + "/" + strconv.Itoa(startlvl) + ".jpg") //открытие изображения как файл
			if err != nil {
				fmt.Println(err.Error())
			}
			img, err := jpeg.Decode(imgfile) //декодирование изображения
			if err != nil {
				fmt.Println(err.Error())
			}
			size := img.Bounds().Size()
			imgfile.Close()
			MatrixPicDataR = TXTtoData("pics/"+picname+"noStr/Daubechi_"+strconv.Itoa(Db)+"/"+strconv.Itoa(startlvl)+".txt", size)
			fmt.Println(Db, "Добеши уровень")
			for k := 0; k < startlvl; k++ {
				cur_part_X := MatrixPicDataR.n / int(math.Pow(2, float64(startlvl-k)))
				cur_part_Y := MatrixPicDataR.m / int(math.Pow(2, float64(startlvl-k)))
				fmt.Println("Размер изображения:", cur_part_X, cur_part_Y)
				RccR := New(cur_part_X, cur_part_Y)
				RcdR := New(cur_part_X, cur_part_Y)
				RdcR := New(cur_part_X, cur_part_Y)
				RddR := New(cur_part_X, cur_part_Y)
				for i := 0; i < cur_part_X*2; i++ {
					for j := 0; j < cur_part_Y*2; j++ {
						if i < cur_part_X {
							if j < cur_part_Y {
								RccR.mat[i][j] = MatrixPicDataR.mat[i][j]
							} else {
								RcdR.mat[i][j-cur_part_Y] = 0 //MatrixPicDataR.mat[i][j]
							}
						} else {
							if j < cur_part_Y {
								RdcR.mat[i-cur_part_X][j] = 0 //MatrixPicDataR.mat[i][j]
							} else {
								RddR.mat[i-cur_part_X][j-cur_part_Y] = 0 //MatrixPicDataR.mat[i][j]
							}
						}
					}
				}

				RccR = RccR.TMat()
				RcdR = RcdR.TMat()
				RdcR = RdcR.TMat()
				RddR = RddR.TMat()
				for i := 0; i < RccR.n; i++ {
					RccR.mat[i] = DaubechiRecoverNormalize(append(RccR.mat[i], RcdR.mat[i]...), Db, 1)
					RdcR.mat[i] = DaubechiRecoverNormalize(append(RdcR.mat[i], RddR.mat[i]...), Db, 1)
				}
				RccR.m *= 2
				RdcR.m *= 2
				RccR = RccR.TMat()
				RdcR = RdcR.TMat()
				for i := 0; i < RccR.n; i++ {
					RccR.mat[i] = DaubechiRecoverNormalize(append(RccR.mat[i], RdcR.mat[i]...), Db, 1)
				}
				RccR.m *= 2
				for i := 0; i < RccR.n; i++ {
					for j := 0; j < RccR.m; j++ {
						MatrixPicDataR.mat[i][j] = RccR.mat[i][j]
					}
				}
				NewData := make([][]float64, 3)
				for i := 0; i < 3; i++ {
					NewData[i] = make([]float64, 0)
				}
				for i := 0; i < MatrixPicDataR.n; i++ {
					NewData[0] = append(NewData[0], MatrixPicDataR.mat[i]...)
				}
				NewData[1], NewData[2] = NewData[0], NewData[0]

				os.MkdirAll("pics/Recover/", 0777)
				os.MkdirAll("pics/Recover/"+picname+"/", 0777)
				os.MkdirAll("pics/Recover/"+picname+"/Daubechi"+strconv.Itoa(Db)+"/", 0777)
				ChangePic("pics/Recover/"+picname+"/Daubechi"+strconv.Itoa(Db)+"/"+strconv.Itoa(startlvl-k)+".jpg", NewData, size)
			}

			goto START
		}
		for Level := 1; ; Level++ {
			fmt.Println("Level = ", Level)
			stop := false
			NewReddata := make([]float64, 0)
			NewGreendata := make([]float64, 0)
			NewBluedata := make([]float64, 0)
			if flag == 1 {
				if HaaraLvlStretching(dataPic[0][0:(1)*size.X], Level) == nil {
					stop = true
					break
				}
				for i := 0; i < size.Y; i++ {
					NewReddata = append(NewReddata, toSlice(HaaraLvlStretching(dataPic[0][i*size.X:(i+1)*size.X], Level))...)
					NewGreendata = append(NewGreendata, toSlice(HaaraLvlStretching(dataPic[1][i*size.X:(i+1)*size.X], Level))...)
					NewBluedata = append(NewBluedata, toSlice(HaaraLvlStretching(dataPic[2][i*size.X:(i+1)*size.X], Level))...)
				}
				if stop {
					fmt.Println("Максимальный уровень сжатия - ", Level-1, ".\nИзображение сжато.")
					break
				}
				picdata := make([][]float64, 3)
				for i := 0; i < size.Y; i++ {
					for j := 0; j < size.X; j++ {
						picdata[0] = append(picdata[0], NewReddata[i+j*size.Y])
						picdata[1] = append(picdata[1], NewGreendata[i+j*size.Y])
						picdata[2] = append(picdata[2], NewBluedata[i+j*size.Y])
					}
				}
				NewReddata = nil
				NewGreendata = nil
				NewBluedata = nil
				if (HaaraLvlStretching((picdata[0][0*size.Y : (0+1)*size.Y]), Level)) == nil {
					stop = true
					fmt.Println("Максимальный уровень сжатия - ", Level-1, ".\nИзображение сжато.")
					break
				}
				for i := 0; i < size.X; i++ {

					NewReddata = append(NewReddata, toSlice(HaaraLvlStretching((picdata[0][i*size.Y:(i+1)*size.Y]), Level))...)
					NewGreendata = append(NewGreendata, toSlice(HaaraLvlStretching((picdata[1][i*size.Y:(i+1)*size.Y]), Level))...)
					NewBluedata = append(NewReddata, toSlice(HaaraLvlStretching((picdata[2][i*size.Y:(i+1)*size.Y]), Level))...)
				}
				if stop {
					fmt.Println("Максимальный уровень сжатия - ", Level-1, ".\nИзображение сжато.")
					break
				}
				bufRed := make([]float64, 0)
				bufGreen := make([]float64, 0)
				bufBlue := make([]float64, 0)
				for i := 0; i < size.X; i++ {
					for j := 0; j < size.Y; j++ {
						bufRed = append(bufRed, NewReddata[i+j*size.X])
						bufGreen = append(bufGreen, NewGreendata[i+j*size.X])
						bufBlue = append(bufBlue, NewBluedata[i+j*size.X])
					}
				}
				NewReddata = bufRed
				NewGreendata = bufGreen
				NewBluedata = bufBlue
			} else if flag == 2 {
				check, _ := DaubechiNormalize(dataPic[0][0*size.X:(0+1)*size.X], Db, Level)
				if check == nil {
					stop = true
					break
				}
				for i := 0; i < size.Y; i++ {
					r, _ := DaubechiNormalize(dataPic[0][i*size.X:(i+1)*size.X], Db, Level)
					g, _ := DaubechiNormalize(dataPic[1][i*size.X:(i+1)*size.X], Db, Level)
					b, _ := DaubechiNormalize(dataPic[2][i*size.X:(i+1)*size.X], Db, Level)
					NewReddata = append(NewReddata, r...)
					NewGreendata = append(NewGreendata, g...)
					NewBluedata = append(NewBluedata, b...)
				}
				if stop {
					fmt.Println("Максимальный уровень сжатия - ", Level-1, ".\nИзображение сжато.")
					break
				}
				picdata := make([][]float64, 3)
				for i := 0; i < size.Y; i++ {
					for j := 0; j < size.X; j++ {
						picdata[0] = append(picdata[0], NewReddata[i+j*size.Y])
						picdata[1] = append(picdata[1], NewGreendata[i+j*size.Y])
						picdata[2] = append(picdata[2], NewBluedata[i+j*size.Y])
					}
				}
				NewReddata = nil
				NewGreendata = nil
				NewBluedata = nil
				check, _ = DaubechiNormalize(picdata[0][0*size.Y:(0+1)*size.Y], Db, Level)
				if check == nil {
					fmt.Println("Максимальный уровень сжатия - ", Level-1, ".\nИзображение сжато.")
					stop = true
					break
				}
				for i := 0; i < size.X; i++ {
					r, _ := DaubechiNormalize(picdata[0][i*size.Y:(i+1)*size.Y], Db, Level)
					g, _ := DaubechiNormalize(picdata[1][i*size.Y:(i+1)*size.Y], Db, Level)
					b, _ := DaubechiNormalize(picdata[2][i*size.Y:(i+1)*size.Y], Db, Level)
					NewReddata = append(NewReddata, r...)
					NewGreendata = append(NewGreendata, g...)
					NewBluedata = append(NewReddata, b...)
				}
				if stop {
					fmt.Println("Максимальный уровень сжатия - ", Level-1, ".\nИзображение сжато.")
					break
				}
				bufRed := make([]float64, 0)
				bufGreen := make([]float64, 0)
				bufBlue := make([]float64, 0)
				for i := 0; i < size.X; i++ {
					for j := 0; j < size.Y; j++ {
						bufRed = append(bufRed, NewReddata[i+j*size.X])
						bufGreen = append(bufGreen, NewGreendata[i+j*size.X])
						bufBlue = append(bufBlue, NewBluedata[i+j*size.X])
					}
				}
				NewReddata = bufRed
				NewGreendata = bufGreen
				NewBluedata = bufBlue
			} else if flag == 3 {
				check, _ := DaubechiNormalizeNoStretching(DataMatrixRed.mat[0], 2, 1)
				if check == nil {
					stop = true
					break
				}
				m_size := DataMatrixRed.m / 2
				DmatR := New(DataMatrixRed.n, m_size)
				xn_size := DataMatrixRed.n
				fmt.Println(DataMatrixRed.n, " ", DataMatrixRed.m)
				for i := 0; i < xn_size; i++ {
					Cr, Dr := DaubechiNormalizeNoStretching(DataMatrixRed.mat[i], 2, 1)
					DataMatrixRed.mat[i] = Cr
					take := len(Dr) / 2
					DmatR.mat[i] = Dr[take:]
				}
				DataMatrixRed.m /= 2
				//Превратили матрицу из NxN в NX(N/2)
				DataMatrixRed = DataMatrixRed.TMat()
				DmatR = DmatR.TMat()
				//Транспонированы в вид (N/2)xN
				m_size = DataMatrixRed.m / 2
				CDr := New(DataMatrixRed.n, m_size)
				DDr := New(DataMatrixRed.n, m_size)
				xn_size = DataMatrixRed.n
				check, _ = DaubechiNormalizeNoStretching(DataMatrixRed.mat[0], 2, 1)
				if check == nil {
					stop = true
					break
				}
				for i := 0; i < xn_size; i++ {
					Cr, Dr := DaubechiNormalizeNoStretching(DataMatrixRed.mat[i], 2, 1)
					DataMatrixRed.mat[i] = Cr //Матрица СС
					take := len(Dr) / 2
					CDr.mat[i] = Dr[take:] //Матрица СD
					DCr1, DDr1 := DaubechiNormalizeNoStretching(DmatR.mat[i], 2, 1)
					DmatR.mat[i] = DCr1      //Матрица DC
					DDr.mat[i] = DDr1[take:] //Матрица DD
				}
				DataMatrixRed.m /= 2
				DmatR.m /= 2
				DataMatrixRed = DataMatrixRed.TMat()
				DmatR = DmatR.TMat()
				DDr = DDr.TMat()
				CDr = CDr.TMat()
				picdara := make([][]float64, 3)
				for i := 0; i < 3; i++ {
					picdara[i] = make([]float64, 0)
				}
				for i := 0; i < DataMatrixRed.n*2; i++ {
					for j := 0; j < DataMatrixRed.m*2; j++ {
						if i < DataMatrixRed.n {
							if j < DataMatrixRed.m {
								MatrixPicDataR.mat[i][j] = DataMatrixRed.mat[i][j]
							} else {
								MatrixPicDataR.mat[i][j] = CDr.mat[i][j-DataMatrixRed.m]
							}
						} else {
							if j < DataMatrixRed.m {
								MatrixPicDataR.mat[i][j] = DmatR.mat[i-DataMatrixRed.n][j]
							} else {
								MatrixPicDataR.mat[i][j] = DDr.mat[i-DataMatrixRed.n][j-DataMatrixRed.m]
							}
						}
					}
				}
				fmt.Println("Matrix size ", len(MatrixPicDataR.mat), " ", len(MatrixPicDataR.mat[0]))
				for i := 0; i < MatrixPicDataR.n; i++ {
					picdara[0] = append(picdara[0], MatrixPicDataR.mat[i]...)
				}
				picdara[1], picdara[2] = picdara[0], picdara[0]
				ChangePic("pics/TestALL_"+strconv.Itoa(Level)+".jpg", picdara, size)

				fmt.Println(DataMatrixRed.n, " ", DataMatrixRed.m)
			} else if flag == 4 {
				check, _ := DaubechiNormalizeNoStretching(DataMatrixRed.mat[0], Db, 1)
				if check == nil {
					stop = true
					break
				}
				m_size := DataMatrixRed.m / 2
				DmatR := New(DataMatrixRed.n, m_size)
				xn_size := DataMatrixRed.n
				fmt.Println(DataMatrixRed.n, " ", DataMatrixRed.m)
				for i := 0; i < xn_size; i++ {
					Cr, Dr := DaubechiNormalizeNoStretching(DataMatrixRed.mat[i], Db, 1)
					DataMatrixRed.mat[i] = Cr
					take := len(Dr) / 2
					DmatR.mat[i] = Dr[take:]
				}
				DataMatrixRed.m /= 2
				//Превратили матрицу из NxN в NX(N/2)
				DataMatrixRed = DataMatrixRed.TMat()
				DmatR = DmatR.TMat()
				//Транспонированы в вид (N/2)xN
				m_size = DataMatrixRed.m / 2
				CDr := New(DataMatrixRed.n, m_size)
				DDr := New(DataMatrixRed.n, m_size)
				xn_size = DataMatrixRed.n
				check, _ = DaubechiNormalizeNoStretching(DataMatrixRed.mat[0], Db, 1)
				if check == nil {
					stop = true
					break
				}
				for i := 0; i < xn_size; i++ {
					Cr, Dr := DaubechiNormalizeNoStretching(DataMatrixRed.mat[i], Db, 1)
					DataMatrixRed.mat[i] = Cr //Матрица СС
					take := len(Dr) / 2
					CDr.mat[i] = Dr[take:] //Матрица СD
					DCr1, DDr1 := DaubechiNormalizeNoStretching(DmatR.mat[i], Db, 1)
					DmatR.mat[i] = DCr1      //Матрица DC
					DDr.mat[i] = DDr1[take:] //Матрица DD
				}
				DataMatrixRed.m /= 2
				DmatR.m /= 2
				DataMatrixRed = DataMatrixRed.TMat()
				DmatR = DmatR.TMat()
				DDr = DDr.TMat()
				CDr = CDr.TMat()
				picdara := make([][]float64, 3)
				for i := 0; i < 3; i++ {
					picdara[i] = make([]float64, 0)
				}
				for i := 0; i < DataMatrixRed.n*2; i++ {
					for j := 0; j < DataMatrixRed.m*2; j++ {
						if i < DataMatrixRed.n {
							if j < DataMatrixRed.m {
								MatrixPicDataR.mat[i][j] = DataMatrixRed.mat[i][j]
							} else {
								MatrixPicDataR.mat[i][j] = CDr.mat[i][j-DataMatrixRed.m]
							}
						} else {
							if j < DataMatrixRed.m {
								MatrixPicDataR.mat[i][j] = DmatR.mat[i-DataMatrixRed.n][j]
							} else {
								MatrixPicDataR.mat[i][j] = DDr.mat[i-DataMatrixRed.n][j-DataMatrixRed.m]
							}
						}
					}
				}
				fmt.Println("Matrix size ", len(MatrixPicDataR.mat), " ", len(MatrixPicDataR.mat[0]))

				for i := 0; i < MatrixPicDataR.n; i++ {
					picdara[0] = append(picdara[0], MatrixPicDataR.mat[i]...)
				}
				picdara[1], picdara[2] = picdara[0], picdara[0]
				srize := size
				srize.X = len(DataMatrixRed.mat)
				srize.Y = len(DataMatrixRed.mat[0])
				picdara1 := make([][]float64, 3)
				for i := 0; i < 3; i++ {
					picdara1[i] = make([]float64, 0)
				}
				for i := 0; i < DataMatrixRed.n; i++ {
					picdara1[0] = append(picdara1[0], DataMatrixRed.mat[i]...)
				}
				picdara1[1], picdara1[2] = picdara1[0], picdara1[0]
				ChangePic("pics/"+picname+"noStr"+"/"+method+"/Pic"+strconv.Itoa(Level)+".jpg", picdara1, srize)
				DataToTXT("pics/"+picname+"noStr"+"/"+method+"/"+strconv.Itoa(Level)+".txt", MatrixPicDataR)
				ChangePic("pics/"+picname+"noStr"+"/"+method+"/"+strconv.Itoa(Level)+".jpg", picdara, size)
				fmt.Println(DataMatrixRed.n, " ", DataMatrixRed.m)
			}
			if stop {
				fmt.Println("Максимальный уровень сжатия - ", Level-1, ".\nИзображение сжато.")
				break
			}
			if (flag == 1) || (flag == 2) {
				os.MkdirAll("pics/"+picname, 0777)
				os.MkdirAll("pics/"+picname+"/"+method, 0777)
				if (NewReddata != nil) && (NewGreendata != nil) && (NewBluedata != nil) {
					NewData := make([][]float64, 3)
					for i := 0; i < 3; i++ {
						NewData[i] = make([]float64, 0)
					}
					NewData[0] = append(NewData[0], (NewReddata)...)
					NewData[1] = append(NewData[1], (NewGreendata)...)
					NewData[2] = append(NewData[2], (NewBluedata)...)
					ChangePic("pics/"+picname+"/"+method+"/"+strconv.Itoa(Level)+".jpg", NewData, size)
				}
			}
		}
		goto START
	case 0:
	default:
		fmt.Println("Неверная выбранная команда")
		goto START
	}
}

func TXTtoData(path string, size image.Point) *Martix {
	file, err := os.Open(path) // создаем файл
	if err != nil {            // если возникла ошибка
		fmt.Println("Unable to create file:", err)
		os.Exit(1) // выходим из программы
	}
	defer file.Close() // закрываем файл
	var content string
	buf := make([]byte, 1024)
	for {
		n, err := file.Read(buf)
		if err == io.EOF { // если конец файла
			break // выходим из цикла
		}
		content += string(buf[:n])
	}
	tokens := strings.FieldsFunc(content, func(r rune) bool {
		return r == ' ' || r == '\n' || r == '\r'
	})
	countY := 0
	countX := 0
	result := New(size.X, size.Y)
	for i := 0; i < len(tokens); i++ {
		number, err := strconv.ParseFloat(tokens[i], 64)
		if err != nil {
			log.Fatal(err)
		}
		result.mat[countX][countY] = number
		countY++
		if countY == size.Y {
			countY = 0
			countX++
		}
	}
	return result
}
func TXTtoData1(path string, size image.Point) *Martix {
	file, err := os.Open(path) // создаем файл
	if err != nil {            // если возникла ошибка
		fmt.Println("Unable to create file:", err)
		os.Exit(1) // выходим из программы
	}
	defer file.Close() // закрываем файл
	var content string
	buf := make([]byte, 1024)
	for {
		n, err := file.Read(buf)
		if err == io.EOF { // если конец файла
			break // выходим из цикла
		}
		content += string(buf[:n])
	}
	var str string
	countY := 0
	countX := 0
	result := New(size.X, size.Y)
	for i := 0; i < len(content); i++ {
		value := (content[i])
		if value == ' ' || value == '\n' || value == '\r' || value == 0 {
			num, err := strconv.ParseFloat(str, 64)
			if err != nil {
				log.Fatal(err)
			}
			result.mat[countX][countY] = num
			countY++
			if countY == size.Y {
				countY = 0
				countX++
			}
			str = ""
		} else {
			str += string(value)
		}
	}
	return result
}

func DataToTXT(path string, data *Martix) {
	file, err := os.Create(path) // создаем файл
	if err != nil {              // если возникла ошибка
		fmt.Println("Unable to create file:", err)
		os.Exit(1) // выходим из программы
	}
	defer file.Close() // закрываем файл
	for i := 0; i < data.n; i++ {
		for j := 0; j < data.m; j++ {
			fmt.Fprint(file, data.mat[i][j])
			fmt.Fprint(file, " ")
		}
		fmt.Fprint(file, "\n")
	}
}
func Equal(first, second *Martix) {
	first.mat = second.mat
	first.m = second.m
	first.n = second.n
}
func Reverse(Red, Green, Blue []float64, newX, newY int) [][]float64 {
	picdata := make([][]float64, 3)
	for i := 0; i < newX; i++ {
		for j := 0; j < newY; j++ {
			picdata[0] = append(picdata[0], Red[i+j*newX])
			picdata[1] = append(picdata[1], Green[i+j*newX])
			picdata[2] = append(picdata[2], Blue[i+j*newX])
		}
	}
	return picdata
}
func DaubechiNormalize(signal []float64, DLevel int, Level int) (res, rec []float64) {
	for i := 0; i < Level; i++ {
		if len(signal)/int(math.Pow(2, float64(i))) < DLevel || len(signal)/int(math.Pow(2, float64(i)))%2 == 1 {
			return nil, nil
		}
	}
	size := len(signal)
	h_k := Dk[DLevel]
	g_k := Generate_gk(h_k)
	res = make([]float64, 0)
	rec = make([]float64, 0)
	buf := make([]float64, 0)
	data := signal
	for k := 0; k < Level; k++ {
		res = nil
		matrix := New(size, 1)
		for i := 0; i < size; i++ {
			matrix.Update(i, 0, data[i])
		}
		arr := Gen_Matrix(h_k, g_k, size).Mul(matrix)
		for i := 0; i < size/2; i++ {
			res = append(res, arr.mat[i][0]/2.0)
		}
		for i := size / 2; i < size; i++ {
			buf = append(buf, arr.mat[i][0]/2.0)
		}
		data = (res)
		rec = append(buf, rec...)
		buf = nil
		if k+1 != Level {
			res = nil
		} else {
			break
		}
		size = len(data)
	}
	rec = append(res, rec...)
	return Stretching(res, Level), rec
}
func DaubechiRecoverNormalize(signal []float64, DLevel int, Level int) []float64 {
	h_k := Dk[DLevel]
	g_k := Generate_gk(h_k)
	size := len(signal)
	data := signal
	for k := 0; k < Level; k++ {
		take := size / int(math.Pow(2, float64(Level-k)))
		matrix := New(take*2, 1)
		matrix.UpdateArray(data[:take*2])
		data = append(data[take*2:], []float64{}...)
		arr := Gen_Matrix(h_k, g_k, take*2).TMat().Mul(matrix).TMat()
		data = append(arr.mat[0], data...)
	}
	return data
}

func Stretching(signal []float64, lvl int) []float64 {
	res := make([]float64, 0)
	for _, v := range signal {
		for i := 0; i < int(math.Pow(2, float64(lvl))); i++ {
			res = append(res, v)
		}
	}
	return res
}
func DaubechiNoStretching(signal []float64, DLevel int, Level int) (res, rec []float64) {
	for i := 0; i < Level; i++ {
		if len(signal)/int(math.Pow(2, float64(i))) < DLevel || len(signal)/int(math.Pow(2, float64(i)))%2 == 1 {
			return nil, nil
		}
	}
	size := len(signal)
	h_k := Dk[DLevel]
	g_k := Generate_gk(h_k)
	res = make([]float64, 0)
	rec = make([]float64, 0)
	buf := make([]float64, 0)
	data := signal
	for k := 0; k < Level; k++ {
		res = nil
		matrix := New(size, 1)
		for i := 0; i < size; i++ {
			matrix.Update(i, 0, data[i])
		}
		arr := Gen_Matrix(h_k, g_k, size)
		arr = arr.Mul(matrix)
		for i := 0; i < size/2; i++ {
			res = append(res, arr.mat[i][0]/math.Sqrt(2))
		}
		for i := size / 2; i < size; i++ {
			buf = append(buf, arr.mat[i][0]/math.Sqrt(2))
		}
		data = (res)
		rec = append(buf, rec...)
		buf = nil
		if k+1 != Level {
			res = nil
		} else {
			break
		}
		size = len(data)
	}
	rec = append(res, rec...)
	return res, rec
}
func DaubechiNormalizeNoStretching(signal []float64, DLevel int, Level int) (res, rec []float64) {

	for i := 0; i < Level; i++ {
		if len(signal)/int(math.Pow(2, float64(i))) < DLevel || len(signal)/int(math.Pow(2, float64(i)))%2 == 1 {
			return nil, nil
		}
	}
	size := len(signal)
	h_k := Dk[DLevel]
	g_k := Generate_gk(h_k)
	res = make([]float64, 0)
	rec = make([]float64, 0)
	buf := make([]float64, 0)
	data := signal
	for k := 0; k < Level; k++ {
		res = nil
		matrix := New(size, 1)
		for i := 0; i < size; i++ {
			matrix.Update(i, 0, data[i])
		}
		arr := Gen_Matrix(h_k, g_k, size).Mul(matrix)
		for i := 0; i < size/2; i++ {
			res = append(res, arr.mat[i][0]/2.0)
		}
		for i := size / 2; i < size; i++ {
			buf = append(buf, arr.mat[i][0]/2.0)
		}
		data = (res)
		rec = append(buf, rec...)
		buf = nil
		if k+1 != Level {
			res = nil
		} else {
			break
		}
		size = len(data)
	}
	rec = append(res, rec...)
	return res, rec
}
func Gen_Matrix(h_k, g_k []float64, size int) *Martix {
	matrix := New(size, size)
	for i := 0; i < size/2; i++ {
		for k := 0; k < len(h_k); k++ {
			if (2*i + k) < size {
				matrix.mat[i][2*i+k] = h_k[k]
				matrix.mat[size/2+i][2*i+k] = g_k[k]
			} else {
				matrix.mat[i][2*i+k-size] = h_k[k]
				matrix.mat[size/2+i][2*i+k-size] = g_k[k]
			}
		}
	}
	return matrix
}

func Generate_gk(h_k []float64) []float64 { //Создание низкочастотных фильтров
	size := len(h_k)
	g_k := make([]float64, size)
	for i := 0; i < size; i++ {
		g_k[i] = math.Pow(-1.0, float64(i)) * h_k[size-i-1]
	}
	return g_k
}

func ChangePic(str string, data [][]float64, size image.Point) {
	m := image.NewRGBA(image.Rect(0, 0, size.X, size.Y))
	counter := 0
	for x := 0; x < size.X; x++ {
		for y := 0; y < size.Y; y++ {
			color := color.RGBA{uint8(data[0][counter]), uint8(data[1][counter]), uint8(data[2][counter]), 255}
			m.Set(x, y, color)
			counter++
		}
	}
	outFile, err := os.Create(str)
	if err != nil {
		log.Fatal(err)
	}
	defer outFile.Close()
	jpeg.Encode(outFile, m, &jpeg.Options{Quality: 100})
}
func GetDataPic_1(str string) ([][]float64, image.Point) {
	imgfile, err := os.Open(str) //открытие изображения как файл
	if err != nil {
		fmt.Println(err.Error())
	}
	defer imgfile.Close()
	img, err := jpeg.Decode(imgfile) //декодирование изображения
	if err != nil {
		fmt.Println(err.Error())
	}
	size := img.Bounds().Size()
	img.ColorModel()
	count := 0
	img.Bounds().Size() //определение границы
	array := make([][]float64, 3)
	for i := 0; i < 3; i++ {
		array[i] = make([]float64, size.X*size.Y)
	}

	for x := 0; x < size.X; x++ {
		for y := 0; y < size.Y; y++ {
			r, g, b, _ := img.At(x, y).RGBA()
			array[0][count] = float64(uint8(r))
			array[1][count] = float64(uint8(g))
			array[2][count] = float64(uint8(b))
			count++
		}
	}
	return array, size
}
func GetDataPic(str string) ([][]float64, image.Point) {
	imgfile, err := os.Open(str) //открытие изображения как файл
	if err != nil {
		fmt.Println(err.Error())
	}
	defer imgfile.Close()
	img, err := jpeg.Decode(imgfile) //декодирование изображения
	if err != nil {
		fmt.Println(err.Error())
	}
	size := img.Bounds().Size()
	count := 0
	img.Bounds().Size() //определение границы
	array := make([][]float64, 3)
	for i := 0; i < 3; i++ {
		array[i] = make([]float64, size.X*size.Y)
	}

	for x := 0; x < size.X; x++ {
		for y := 0; y < size.Y; y++ {
			r, g, b, _ := img.At(x, y).RGBA()
			array[0][count] = float64(uint8(r))
			array[1][count] = float64(uint8(g))
			array[2][count] = float64(uint8(b))
			count++
		}
	}
	return array, size
}

func AddToPlot(Xn []float64, Yn []float64, p *plot.Plot, color color.RGBA) { // добавляет на график точки
	Points := Points(Xn, Yn)
	filled, err := plotter.NewLine(Points)
	filled.Color = color
	if err != nil {
		log.Panic(err)
	}
	p.Add(filled)
	err = p.Save(1000, 600, Path)
	if err != nil {
		log.Panic(err)
	}
}

func CreatePlot(title string, w, h font.Length) *plot.Plot {
	p := plot.New()
	p.Title.Text = title
	//p.X.Label.Text = "X"
	//p.Y.Label.Text = "Y"

	p.Title.TextStyle.Font.Size = vg.Points(32) // Размер заголовка
	p.X.Tick.Label.Font.Size = vg.Points(32)    // Размер меток по оси X
	p.Y.Tick.Label.Font.Size = vg.Points(32)
	p.Add(plotter.NewGrid())
	err := p.Save(w, h, Path)
	if err != nil {
		log.Panic(err)
	}
	return p
}

func Points(Xn []float64, Yn []float64) plotter.XYs {
	pts := make(plotter.XYs, len(Xn))
	for i := range pts {
		pts[i].X = Xn[i]
		pts[i].Y = Yn[i]
	}
	return pts
}
func HaaraLvlStretching(Xn []float64, level int) [][]float64 { // Основная функция
	for i := 0; i < level; i++ {
		if len(Xn)/int(math.Pow(2, float64(i)))%2 == 1 && i != level {
			return nil
		}
	}
	size := math.Pow(2, float64(level))
	var array [][]float64 = make([][]float64, len(Xn))
	ksi := 0
	j := 0
	for i := 0; i < len(array); i++ {
		array[i] = make([]float64, 0)
		array[i] = append(array[i], Haara(Xn[ksi*(int(size)):ksi*(int(size))+(int(size))], level)...)
		for j = 1; j < int(size); j++ {
			array[i+j] = make([]float64, 1)
			array[i+j] = array[i]
		}
		ksi++
		i += j - 1
	}
	return array
}
func HaaraTr(Xn []float64, level int) []float64 { //обратный подсчет коэфицентов вейвлет преобразования
	sizeSignal := len(Xn)
	if sizeSignal&(sizeSignal-1) != 0 {
		fmt.Println("Сигнал не кратен 2^N!")
		return nil
	}
	HaarMatrix := New(sizeSignal, sizeSignal)
	for i := 0; i < sizeSignal; i++ {
		HaarMatrix.Update(0, i, 1.0)
		if i < sizeSignal/2 {
			HaarMatrix.Update(1, i, 1.0)
		} else {
			HaarMatrix.Update(1, i, -1.0)
		}
	}
	i := 2
	for step := 1; step < int(math.Log2(float64(sizeSignal))); step++ {
		part := math.Pow(2, float64(step))
		startCol := 0
		for z := 0; z < int(part); z++ {
			for j := 0; j < sizeSignal/int(part); j++ {
				HaarMatrix.Update(i+z, startCol+j, math.Sqrt(math.Pow(2, float64(step)))*MomWavelet(HaarMatrix.mat[1], int(math.Pow(2, float64(step)))*(startCol+j)-sizeSignal*z))
			}
			startCol += sizeSignal / int(part)
		}
		i += int(part)
	}
	HaarMatrix.Mulconst(1.0 / math.Sqrt(float64(sizeSignal)))
	HaarMatrix.Mulconst(math.Sqrt(math.Pow(2.0, float64(level))))
	HaarMatrix.Mulconst(1.0 / float64(sizeSignal))
	HaarMatrix.Mulconst(math.Sqrt(float64(sizeSignal)))
	HaarMatrix = HaarMatrix.TMat()
	vector := New(len(Xn), 1)
	vector.UpdateArray(Xn)
	return toSlice(HaarMatrix.Mul(vector).mat)
}
func HaaraLvlStretchingNormalize(Xn []float64, level int) [][]float64 { // Основная функция
	for i := 0; i < level; i++ {
		if len(Xn)/int(math.Pow(2, float64(i)))%2 == 1 && i != level {
			fmt.Println(math.Pow(2, float64(level)))
			return nil
		}
	}
	size := math.Pow(2, float64(level))
	var array [][]float64 = make([][]float64, len(Xn))
	ksi := 0
	j := 0
	for i := 0; i < len(array); i++ {
		array[i] = make([]float64, 0)
		array[i] = append(array[i], HaaraNormalize(Xn[ksi*(int(size)):ksi*(int(size))+(int(size))], level)...)
		for j = 1; j < int(size); j++ {
			array[i+j] = make([]float64, 1)
			array[i+j] = array[i]
		}
		ksi++
		i += j - 1
	}
	return array
}
func HaaraNormalize(Xn []float64, level int) []float64 { //вейвлет преобразование
	sizeSignal := len(Xn)
	if len(Xn)/int(math.Pow(2, float64(level)))%2 != 1 {
		fmt.Println("Сигнал не кратен 2^N")
		return nil
	}
	HaarMatrix := New(sizeSignal, sizeSignal)
	for i := 0; i < sizeSignal; i++ {
		HaarMatrix.Update(0, i, 1.0)
		if i < sizeSignal/2 {
			HaarMatrix.Update(1, i, 1.0)
		} else {
			HaarMatrix.Update(1, i, -1.0)
		}
	}
	i := 2
	for step := 1; step < int(math.Log2(float64(sizeSignal))); step++ {
		part := math.Pow(2, float64(step))
		startCol := 0
		for z := 0; z < int(part); z++ {
			for j := 0; j < sizeSignal/int(part); j++ {
				HaarMatrix.Update(i+z, startCol+j, math.Sqrt(math.Pow(2, float64(step)))*MomWavelet(HaarMatrix.mat[1], int(math.Pow(2, float64(step)))*(startCol+j)-sizeSignal*z))
			}
			startCol += sizeSignal / int(part)
		}
		i += int(part)
	}
	HaarMatrix.Mulconst(1.0 / math.Sqrt(float64(sizeSignal)))
	HaarMatrix.Mulconst(1.0 / math.Sqrt(math.Pow(2.0, float64(level))))
	vector := New(len(Xn), 1)
	vector.UpdateArray(Xn)
	str := toSlice(HaarMatrix.Mul(vector).mat)
	return str
}
func Haara(Xn []float64, level int) []float64 { //вейвлет преобразование
	sizeSignal := len(Xn)
	if len(Xn)/int(math.Pow(2, float64(level)))%2 != 1 {
		fmt.Println("Сигнал не кратен 2^N")
		return nil
	}
	HaarMatrix := New(sizeSignal, sizeSignal)
	for i := 0; i < sizeSignal; i++ {
		HaarMatrix.Update(0, i, 1.0)
		if i < sizeSignal/2 {
			HaarMatrix.Update(1, i, 1.0)
		} else {
			HaarMatrix.Update(1, i, -1.0)
		}
	}
	i := 2
	for step := 1; step < int(math.Log2(float64(sizeSignal))); step++ {
		part := math.Pow(2, float64(step))
		startCol := 0
		for z := 0; z < int(part); z++ {
			for j := 0; j < sizeSignal/int(part); j++ {
				HaarMatrix.Update(i+z, startCol+j, math.Sqrt(math.Pow(2, float64(step)))*MomWavelet(HaarMatrix.mat[1], int(math.Pow(2, float64(step)))*(startCol+j)-sizeSignal*z))
			}
			startCol += sizeSignal / int(part)
		}
		i += int(part)
	}
	HaarMatrix.Mulconst(1.0 / math.Sqrt(float64(sizeSignal)))
	HaarMatrix.Mulconst(1.0 / math.Sqrt(math.Pow(2.0, float64(level))))
	HaarMatrix.Mulconst(float64(sizeSignal))
	HaarMatrix.Mulconst(1.0 / math.Sqrt(float64(sizeSignal)))
	vector := New(len(Xn), 1)
	vector.UpdateArray(Xn)
	str := toSlice(HaarMatrix.Mul(vector).mat)
	return str
}

func toSlice(mat [][]float64) []float64 {
	slice := make([]float64, 0)
	for _, v := range mat {
		slice = append(slice, v[0])
	}
	return slice
}

func StartCol(mat Martix, row int) int {
	count := 0
	for j := 0; j < mat.m; j++ {
		if mat.mat[row][j] != 0.0 {
			count++
		} else {
			return count
		}
	}
	return 0
}

func MomWavelet(wave []float64, k int) float64 {
	if (k < 0) || (k >= len(wave)) {
		return 0
	} else {
		return wave[k]
	}
}

func OpenFile(name string, N int) (*[]coord, error) { // name - имя файла , N - кол-во символов
	file, err := os.Open(name) // открытие файла
	if err != nil {
		log.Fatal(err)
		os.Exit(1)
	}
	defer file.Close()
	//Чтение данных
	var arr []coord = ReadData(file)
	arr = arr[:N]
	return &arr, nil
}

func ReadData(f *os.File) (arr []coord) {
	var x, y float64
	for i := 0; ; i++ {
		_, err := fmt.Fscanf(f, "%f %f\n", &x, &y)
		arr = append(arr, coord{x: x, y: y})
		if err != nil {
			if err == io.EOF {
				break
			} else {
				fmt.Println(err)
				os.Exit(1)
			}
		}

	}
	f.Seek(0, 0)
	return arr
}

func SignalCreate(size, parts float64) string {
	name := "Signal.txt"
	file, err := os.Create(name)
	if err != nil {
		log.Fatal(err)
		os.Exit(1)
	}
	defer file.Close()
	for x := 0.0; x < size; x += size / parts {
		if x < size/4 {
			fmt.Fprintf(file, "%-.2f %-.5f\n", x, math.Cos(4*pi*10*x))
		} else if x < size/2 {
			fmt.Fprintf(file, "%-.2f %-.5f\n", x, math.Cos(4*pi*25*x))
		} else if x < 3*size/4 {
			fmt.Fprintf(file, "%-.2f %-.5f\n", x, math.Cos(4*pi*50*x))
		} else {
			fmt.Fprintf(file, "%-.2f %-.5f\n", x, math.Cos(2*2*pi*100*x))
		}
	}
	return name
}
func SignalCreateStat(size, parts float64) string {
	name := "SignalStat.txt"
	file, err := os.Create(name)
	if err != nil {
		log.Fatal(err)
		os.Exit(1)
	}
	defer file.Close()
	for x := 0.0; x <= size; x += size / parts {
		fmt.Fprintf(file, "%-.2f %-.5f\n", x, math.Cos(2*pi*10*x)+math.Cos(2*pi*25*x)+math.Cos(2*pi*50*x)+math.Cos(2*pi*100*x))
	}
	return name
}
