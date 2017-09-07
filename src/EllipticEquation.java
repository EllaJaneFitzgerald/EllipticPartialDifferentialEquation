
public class EllipticEquation {
    private final int N;                // задает сетку на области
    private final double h;             // шаг
    private final double eps;           // задает точность приближения
    private final double [][] uReal;    // точное решение

    public EllipticEquation(int m, double sh, double e){
        N = m;
        h = sh;
        eps = e;
        uReal = new double[N][N];

        for(int i=0; i<N; i++){
            for (int j=0; j<N; j++){
                uReal[i][j] = fi(i*h,j*h);
            }
        }
    }

    // задача вида
    // -d/dx1(K1(x1,x2) du/dx1) - d/dx2(K2(x1,x2) du/dx2) = f(x1,x2), (x1,x2) in G
    // u = fi(x1,x2), (x1,x2) in B(G)
    private double k1(double x1){
        return(1);
    }
    private double k2(double x2){
        double r= 1+Math.sqrt(1+x2*x2);
        return(r);
    }
    private double fi(double x1, double x2){            // условие на границе куба
        return(x1*x1*x1+(x1*x1+1)*x2);
    }
    private double f(double x1,double x2){
        return(-6*x1-2*x2-x2*(x1*x1+1)/Math.sqrt(x2*x2+1));
    }

    public boolean Check(double[][] u){                  // проверяет, достаточно ли приближение
        boolean flag = true;
        for (int i=0; i<N; i++){
            for(int j=0;j<N;j++){
                if (Math.abs(u[i][j]-uReal[i][j])>=eps){flag=false;}
            }
        }
        return(!flag);
    }

    public void write(double[][] u){                     // выводит на экран матрицу u
        for (int i=0; i<N;i++){
            for (int j=0;j<N; j++){
                System.out.printf("%.6f",u[i][j]);
                System.out.print("  ");
            }
            System.out.println();
        }
        System.out.println();
    }

    public void methodSimpleIteration(){                                  // метод простой итерации
        double lambdaMin = 8/(h*h)*Math.sin(Math.PI*h/2)*Math.sin(Math.PI*h/2);                                 // минимальное собственное число
        double lambdaMax = (4/(h*h) + (1+Math.sqrt(2))*4/(h*h))*Math.cos(Math.PI*h/2)*Math.cos(Math.PI*h/2);    // максимальное собственное число
        double tau = 2/(lambdaMin+lambdaMax);

        int col = 0;
        double[][] u = new double[N][N];
        double[][] u2 = new double[N][N];
        for(int i=0; i<N; i++){
            for (int j=0; j<N;j++){
                if ((i==0)||(i==N-1)||(j==0)||(j==N-1)){
                    u[i][j] = fi(i*h,j*h);
                    u2[i][j] = fi(i*h,j*h);
                } else {
                    u[i][j] = 1;
                    u2[i][j] = 1;
                }
            }
        }
        while (Check(u)){
            for (int i=1; i<(N-1); i++){
                for(int j=1; j<(N-1);j++){
                    u2[i][j] = +u[i-1][j]*k1(h*(i-0.5))+u[i][j-1]*k2(h*(j-0.5));
                    u2[i][j] = u2[i][j]-u[i][j]*(k1(h*(i-0.5))+k1(h*(i+0.5))+k2(h*(j-0.5))+k2(h*(j+0.5)));
                    u2[i][j] = u2[i][j]+u[i][j+1]*k2(h*(j+0.5))+u[i+1][j]*k1(h*(i+0.5));
                    u2[i][j] = (u2[i][j]/(h*h)+f(i*h,j*h))*tau+u[i][j];
                }
            }
            for (int k=0; k<N; k++){
                for (int l=0; l<N; l++){
                    u[k][l] = u2[k][l];
                }
            }
            col++;
        }
        System.out.println("МЕТОД ПРОСТОЙ ИТЕРАЦИИ");
        System.out.println("Приближенное решение");
        write(u);
        System.out.println("Точное решение");
        write(uReal);
        System.out.println("Количество итераций: "+col);
    }

    public void methodNekrasov(){                                         // метод Некрасова
        int col = 0;
        double[][] u = new double[N][N];
        for(int i=0; i<N; i++){
            for (int j=0; j<N;j++){
                if ((i==0)||(i==N-1)||(j==0)||(j==N-1)){
                    u[i][j] = fi(i*h,j*h);
                } else {
                    u[i][j] = 1;
                }
            }
        }
        double s=0;
        while (Check(u)){
            for (int i=1; i<(N-1); i++){
                for(int j=1; j<(N-1);j++){
                    u[i][j] = k1(h*(i-0.5))*u[i-1][j]+k1(h*(i+0.5))*u[i+1][j]+k2(h*(j-0.5))*u[i][j-1]+k2(h*(j+0.5))*u[i][j+1];
                    u[i][j] = u[i][j]/(h*h) + f(i*h,j*h);
                    s = k1(h*(i-0.5))+k1(h*(i+0.5))+k2(h*(j-0.5))+k2(h*(j+0.5));
                    //s = s/(h*h);
                    u[i][j] = u[i][j]*(h*h)/s;
                }
            }
            col++;
        }
        System.out.println("МЕТОД НЕКРАСОВА");
        System.out.println("Приближенное решение");
        write(u);
        System.out.println("Точное решение");
        write(uReal);
        System.out.println("Количество итераций: "+col);
    }

    public void methodRelaxation(){                                                // метод верхней релаксации
        int col = 0;
        double omega = 1.6;                         // параметр
        double[][] u = new double[N][N];
        double[][] u2 = new double[N][N];
        for(int i=0; i<N; i++){
            for (int j=0; j<N;j++){
                if ((i==0)||(i==N-1)||(j==0)||(j==N-1)){
                    u[i][j] = fi(i*h,j*h);
                    u2[i][j] = fi(i*h,j*h);
                } else {
                    u2[i][j] = 1;
                    u[i][j] = 1;
                }
            }
        }
        double s=0;
        while (Check(u)){
            for (int i=1; i<(N-1); i++){
                for(int j=1; j<(N-1);j++){
                    u2[i][j] = k1(h*(i+0.5))*(u[i+1][j]-u[i][j]) - k1(h*(i-0.5))*(u[i][j]-u2[i-1][j]);
                    u2[i][j] = u2[i][j] + k2(h*(j+0.5))*(u[i][j+1]-u[i][j]) - k2(h*(j-0.5))*(u[i][j]-u2[i][j-1]);
                    u2[i][j] = u2[i][j]/(h*h)+f(i*h,j*h);
                    s = k1(h*(i-0.5))+k1(h*(i+0.5))+k2(h*(j-0.5))+k2(h*(j+0.5));
                    s = s/(h*h);
                    u2[i][j] = u2[i][j]/s;
                    u2[i][j] = u[i][j]+omega*u2[i][j];
                }
            }
            for (int k=0; k<N; k++){
                for (int l=0; l<N; l++){
                    u[k][l] = u2[k][l];
                }
            }
            col++;
        }
        System.out.println("МЕТОД ВЕРХНЕЙ РЕЛАКСАЦИИ");
        System.out.println("Приближенное решение");
        write(u);
        System.out.println("Точное решение");
        write(uReal);
        System.out.println("Количество итераций: "+col);
    }
}
