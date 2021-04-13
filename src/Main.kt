import kotlin.math.abs
import kotlin.math.max

fun main() {
    val firsttask: Array<Array<Double>> = Array(4, { Array(4, {0.0}) })
    firsttask[0] = arrayOf(1.0,0.55,-0.13,0.34)
    firsttask[1] = arrayOf(0.13,-0.17,0.33,0.17)
    firsttask[2] = arrayOf(0.11,0.18,-0.22,-0.11)
    firsttask[3] = arrayOf(0.13,-0.12,0.21,0.22)

    val aminusone: Array<Array<Double>> = Array(4, { Array(4, {0.0}) })
    aminusone[0] = arrayOf(-1.53034,4.08854,12.0119,5.2117)
    aminusone[1] = arrayOf(4.45859,-1.68133,-20.0751,-15.6289)
    aminusone[2] = arrayOf(2.32367,4.46726,-11.3646,-12.7254)
    aminusone[3] = arrayOf(1.11821,-7.59725,-7.19997,5.08793)

    val a1minusone: Array<Array<Double>> = Array(4, { Array(4, {0.0}) })
    a1minusone[0] = arrayOf(1.0,0.0,0.0,0.0)
    a1minusone[1] = arrayOf(0.764706,-5.88235,0.0,0.0)
    a1minusone[2] = arrayOf(1.12567,-4.81283,-4.54545,0.0)
    a1minusone[3] = arrayOf(-1.2483,1.38551,4.33884,4.54545)

    val a2: Array<Array<Double>> = Array(4, { Array(4, {0.0}) })
    a2[0] = arrayOf(0.0,0.55,-0.13,0.34)
    a2[1] = arrayOf(0.0,0.0,0.33,0.17)
    a2[2] = arrayOf(0.0,0.0,0.0,-0.11)
    a2[3] = arrayOf(0.0,0.0,0.0,0.0)

    val numbers: Array<Double> = arrayOf(0.13,0.11, 1.0, 0.18)



    val aa: Array<Array<Double>> = Array(4, { Array(4, {0.0}) })

    val frt: Array<Array<Double>> = Array(4, { Array(4, {0.0}) })
    for (i in 0..3){
        for(k in 0..3){
            frt[i][k] = firsttask[k][i]
        }
    }
    println("a*at")
    for (i in 0..3){
        for(k in 0..3){
            aa[i][k] = frt[i][1]*firsttask[1][k]+frt[i][2]*firsttask[2][k]+frt[i][3]*firsttask[3][k]+frt[i][0]*firsttask[0][k]
            print(aa[i][k])
            print(" ")
        }
        println(" ")
    }
    val numbers1: Array<Double> = matrnavec(frt,numbers)
    print("ynew = ")
    for (i in 0..3){
        print(numbers1[i])
        print(" ")
    }
    println(" ")
    print("определитель для 3 теоремы = ")
    println(det4(aa))
    val x0:Array<Double> = arrayOf(0.0,0.0,0.0,0.0)
    var res:Array<Double> = schit(x0,numbers1,aa)
    print("Х1 = ")
    showvect(res)
    print("Третья норма для р0 = ")
    print(thirdnorm(numbers1))
    var x:Double =  0.38507/0.00088
    var cnt:Int = 0
    while (x>0.0001){
        x=x*(1.47348-0.00088)/(1.47348+0.00088)
        cnt++
    }
    println(" ")
    print(cnt)
    for (i in 0..12803){
        res = schit(res,numbers1,aa)
    }
    println(" ")
    print("Х* = ")
    showvect(res)
    cnt = 0
    x = 17.71428
    while (x>0.001){
        x*=0.72
        cnt++
    }
    println(" ")
    print(cnt)
    val A:Array<Array<Double>> = Array(4, { Array(4, {0.0}) })
    val b: Array<Double> = arrayOf(-1.42,0.48, -2.34, 0.72)
    A[0] = arrayOf(0.17,0.27,-0.13,-0.11)
    A[1] = arrayOf(0.13,-0.12,0.09,-0.06)
    A[2] = arrayOf(0.11,0.05,-0.02,0.12)
    A[3] = arrayOf(0.13,0.18,0.24,0.43)
    showvect(mpi(x0,b,A,cnt))
    println("Сейчас будет МПИ!:")
    /*val ay: Array<Double> = arrayOf(13.181,-22.412,-12.892, -6.975)
    val tau:Array<Array<Double>> = Array(4, { Array(4, {0.0}) })
    tau[0] = arrayOf(-0.016,-0.031,0.033,0.006)
    tau[1] = arrayOf(0.036,0.036,-0.024,0.008)
    tau[2] = arrayOf(0.02,0.028,-0.027,-0.007)
    tau[3] = arrayOf(0.017,0.011,-0.004,0.004)

    var xnach: Array<Double> = arrayOf(13.181,-22.412,-12.892, -6.975)
    for (i in 0..5){
        xnach = plusvec(matrnavec(tau,xnach),ay)
        showvect(xnach)
        println()
    }*/



}
fun det4(A: Array<Array<Double>>): Double {
    return A[0][0]*(A[1][1]*A[2][2]*A[3][3]+A[1][2]*A[2][3]*A[3][1]+A[3][2]*A[1][3]*A[2][1]-A[1][3]*A[2][2]*A[3][1]-A[1][1]*A[3][2]*A[2][3]-A[1][2]*A[2][1]*A[3][3])-
            A[1][0]*(A[0][1]*A[2][2]*A[3][3]+A[0][2]*A[2][3]*A[3][1]+A[3][2]*A[0][3]*A[2][1]-A[0][3]*A[2][2]*A[3][1]-A[0][1]*A[3][2]*A[2][3]-A[0][2]*A[2][1]*A[3][3])+
            A[2][0]*(A[0][1]*A[1][2]*A[3][3]+A[0][2]*A[1][3]*A[3][1]+A[3][2]*A[0][3]*A[1][1]-A[0][3]*A[1][2]*A[3][1]-A[0][1]*A[3][2]*A[1][3]-A[0][2]*A[1][1]*A[3][3])-
            A[3][0]*(A[0][1]*A[1][2]*A[2][3]+A[0][2]*A[1][3]*A[2][1]+A[2][2]*A[0][3]*A[1][1]-A[0][3]*A[1][2]*A[2][1]-A[0][1]*A[2][2]*A[1][3]-A[0][2]*A[1][1]*A[2][3])
}
fun firtsnorm(A: Array<Array<Double>>):Double{
    var x1:Double = 0.0
    var x2:Double = 0.0
    var x3:Double = 0.0
    var x4:Double = 0.0
    for (i in 0..3){
        x1 += abs(A[0][i])
    }
    for (i in 0..3){
        x2 += abs(A[1][i])
    }
    for (i in 0..3){
        x3 += abs(A[2][i])
    }
    for (i in 0..3){
        x4 += abs(A[3][i])
    }
    val x = max(max(x1,x2), max(x3,x4))
    return x
}
fun scalar(A:Array<Double>,B:Array<Double>):Double{
    return A[0]*B[0]+ A[1]*B[1]+A[3]*B[3]+A[2]*B[2]
}
fun matrnavec(A: Array<Array<Double>>,B:Array<Double>):Array<Double>{
    val r: Array<Double> = arrayOf(0.0,0.0, 0.0, 0.0)
    for (i in 0..3){
        r[i] = A[i][0]*B[0]+A[i][1]*B[1]+A[i][2]*B[2]+A[i][3]*B[3]
    }
    return r
}
fun plusvec(x:Array<Double>,b:Array<Double>):Array<Double>{
    val r: Array<Double> = arrayOf(0.0,0.0, 0.0, 0.0)
    for (i in 0..3){
        r[i] = x[i]+b[i]
    }
    return r
}
fun scalumn(x:Array<Double>,b:Double):Array<Double>{
    val r: Array<Double> = arrayOf(0.0,0.0, 0.0, 0.0)
    for (i in 0..3){
        r[i] = x[i]*b
    }
    return r
}
fun schit(x:Array<Double>,b:Array<Double>,A: Array<Array<Double>>):Array<Double>{
    val r: Array<Double> = arrayOf(0.0,0.0, 0.0, 0.0)
    val rk: Array<Double> = plusvec(b,matrnavec(A,scalumn(x,(-1.0))))
    val tauk: Double = scalar(rk,rk)/scalar(rk,matrnavec(A,rk))
    for (i in 0..3){
        r[i] = x[i]+tauk*rk[i]
    }


    return r
}
fun showvect(x:Array<Double>){
    for (i in 0..3){
        print(x[i])
        print(" ")
    }
}
fun thirdnorm(x: Array<Double>):Double{
    return Math.sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3])
}
fun mpi(x:Array<Double>,b:Array<Double>,A: Array<Array<Double>>,k:Int):Array<Double>{
    var r: Array<Double> = x
    var cnt = 0

    while (cnt<k){
        r=plusvec(b,matrnavec(A,r))
        cnt++
    }
    return r
}