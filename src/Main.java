public class Main {

    public static void main(String[] args) {
        int N = 11;
        double h = 0.1;
        double eps = 0.0001;
        EllipticEquation el = new EllipticEquation(N,h,eps);

        el.methodSimpleIteration();
        System.out.println();
        System.out.println();
        el.methodNekrasov();
        System.out.println();
        System.out.println();
        el.methodRelaxation();
    }
}
