import java.util.Random;

import Plotter.Plotter;

public class SortingComparisons {

    final static int INSERTION_VS_QUICK_LENGTH = 12;
    final static int MERGE_VS_QUICK_LENGTH = 15;
    final static int INSERTION_VS_QUICK_SORTED_LENGTH = 12;
    final static int ARBITRARY_VS_RANDOM_LENGTH = 16;
    final static int COUNTING_VS_QUICK_LENGTH = 15;
    final static double T = 600.0;

    /**
     * Sorts a given array using the quick sort algorithm. At each stage the pivot
     * is chosen to be the rightmost element of the subarray.
     * 
     * Should run in average complexity of O(nlog(n)), and worst case complexity of
     * O(n^2)
     * 
     * @param arr - the array to be sorted
     */
    public static void quickSortArbitraryPivot(double[] arr) {
        quickSort(arr, 0, arr.length - 1);
    }

    /**
     * Sorts a given array using the quick sort algorithm. At each stage the pivot
     * is chosen to be the rightmost element of the subarray.
     *
     * @param arr   - the array to be sorted
     * @param start - the first index of subarray should be sorted
     * @param end   - the last index of subarray should be sorted
     */
    public static void quickSort(double[] arr, int start, int end) {
        if (start < end) {
            int q = partition(arr, start, end);
            quickSort(arr, start, q - 1);
            quickSort(arr, q + 1, end);
        }
    }

    /**
     * chosen to be the rightmost element of the subarray to be the pivot. Take all
     * the small elements than the pivot to be at the left side of the array and all
     * the big elements than the pivot on the other side.
     * 
     * @param arr   - the array to be sorted
     * @param start - the first index of subarray
     * @param end   - the last index of the subarray
     * @return the index of pivot after the sort
     */
    public static int partition(double[] arr, int start, int end) {
        double pivot = arr[end];
        int i = start - 1;
        for (int j = start; j < end; j++) {
            if (arr[j] <= pivot) {
                i++;
                swap(arr, i, j);
            }
        }
        if (arr[end] < arr[i + 1]) {
            swap(arr, i + 1, end);
        }
        return i + 1;
    }

    /**
     * swap between two elements of the given array
     * 
     * @param arr - the array should be change elements
     * @param i   - the first index of the element should be replace
     * @param j   - the second index of the element should be replace
     */
    public static void swap(double[] arr, int i, int j) {
        double temp = arr[i];
        arr[i] = arr[j];
        arr[j] = temp;
    }

    /**
     * Sorts a given array using the quick sort algorithm. At each stage the pivot
     * is chosen in the following way: Choose a random index from the range, the
     * element at this index is the pivot.
     * 
     * Should run in average complexity of O(nlog(n)), and worst case complexity of
     * O(n^2)
     * 
     * @param arr - the array to be sorted
     */
    public static void quickSortRandomPivot(double[] arr) {
        quickSortRand(arr, 0, arr.length - 1);
    }

    /**
     * Sorts a given array using the quick sort algorithm. At each stage the pivot
     * is chosen to be a random element of the subarray.
     * 
     * @param arr   - the array to be sorted
     * @param start - the first index of subarray should be sorted
     * @param end   - the last index of subarray should be sorted
     */
    public static void quickSortRand(double[] arr, int start, int end) {
        if (start < end) {
            int q = partitionRand(arr, start, end);
            quickSortRand(arr, start, q - 1);
            quickSortRand(arr, q + 1, end);
        }
    }

    /**
     * chosen a random element of the subarray to be the pivot and swap him with the
     * rightmost element. Then call to partition function.
     * 
     * @param arr   - the array to be sorted
     * @param start - the first index of subarray
     * @param end   - the last index of the subarray
     * @return the index of pivot after the sort
     */
    public static int partitionRand(double[] arr, int start, int end) {
        int x = (int) (start + (Math.random() * (end - start + 1)));
        double pivot = arr[x];
        swap(arr, x, end);
        return partition(arr, start, end);
    }

    /**
     * Sorts a given array using the merge sort algorithm.
     * 
     * Should run in complexity O(nlog(n)) in the worst case.
     * 
     * @param arr - the array to be sorted
     */
    public static void mergeSort(double[] arr) {
        mergeSort(arr, 0, arr.length - 1);
    }

    /**
     * Sorts a given array using the merge sort algorithm.
     * 
     * Should run in complexity O(nlog(n)) in the worst case.
     * 
     * @param arr   - the array to be sorted
     * @param start - the first index of the subarray
     * @param end   - the last index of the subarray
     */
    public static void mergeSort(double[] arr, int start, int end) {
        if (start < end) {
            int mid = (start + end) / 2;
            mergeSort(arr, start, mid);
            mergeSort(arr, mid + 1, end);
            merge(arr, start, mid, end);
        }
    }

    /**
     * Divide the array to two sorted subarrays,from start to mid and from mid to
     * the end/ Proudce one array sorted
     * 
     * @param arr   - the array to be sorted
     * @param start - the first index of the subarray
     * @param mid   - the middle index of the subarray
     * @param end   - the last index of the subarray
     */
    public static void merge(double[] arr, int start, int mid, int end) {
        int n1 = mid - start + 1;
        int n2 = end - mid;
        double[] arrLeft = new double[n1 + 1];
        double[] arrRight = new double[n2 + 1];
        for (int i = 0; i < n1; i++) {
            arrLeft[i] = arr[start + i];
        }
        for (int j = 0; j < n2; j++) {
            arrRight[j] = arr[mid + j + 1];
        }
        double x = Double.MAX_VALUE;
        arrLeft[n1] = x;
        arrRight[n2] = x;
        int i = 0;
        int j = 0;
        for (int k = start; k < end + 1; k++) {
            if (arrLeft[i] <= arrRight[j]) {
                arr[k] = arrLeft[i];
                i++;
            } else {
                arr[k] = arrRight[j];
                j++;
            }
        }
    }

    /**
     * Sorts a given array, using the counting sort algorithm. You may assume that
     * all elements in the array are between 0 and k (not including k).
     * 
     * Should run in complexity O(n + k) in the worst case.
     * 
     * @param arr
     * @param k   - an upper bound for all elements in the array.
     */
    public static void countingSort(int[] arr, int k) {
        int[] temp = new int[k];
        for (int i = 0; i < arr.length; i++) {
            temp[arr[i]]++;
        }
        int l = 0;
        for (int i = 0; i < temp.length; i++) {
            int x = temp[i];
            for (int j = 0; j < x; j++) {
                arr[l] = i;
                l++;
            }
        }
    }

    /**
     * Sorts a given array using insertion sort.
     * 
     * The algorithm should run in complexity O(n^2) in the worst case.
     * 
     * @param arr - the array to be sorted
     */
    public static void insertionSort(double[] arr) {
        for (int i = 1; i < arr.length; i++) {
            double key = arr[i];
            int j = i - 1;
            while (j >= 0 && arr[j] > key) {
                arr[j + 1] = arr[j];
                j--;
            }
            arr[j + 1] = key;
        }

    }

    public static void main(String[] args) {
        int[] arr = {10,9,8,7,6,5,4,3,2,18};
        int k = 19;
        countingSort(arr,k);
        for(int i = 0; i < arr.length; i++) {
            System.out.print(arr[i] + ",");
        }
        //insertionVsQuick();
        //mergeVsQuick();
        //insertionVsQuickOnSortedArray();
        //countingVsQuick();
        //arbitraryPivotVsRandomPivot();
    }

    private static void countingVsQuick() {
        double[] quickTimes = new double[COUNTING_VS_QUICK_LENGTH];
        double[] countingTimes = new double[COUNTING_VS_QUICK_LENGTH];
        long startTime, endTime;
        Random r = new Random();
        for (int i = 0; i < COUNTING_VS_QUICK_LENGTH; i++) {
            long sumQuick = 0;
            long sumCounting = 0;
            for (int k = 0; k < T; k++) {
                int size = (int) Math.pow(2, i);
                double[] a = new double[size];
                int[] b = new int[size];
                for (int j = 0; j < a.length; j++) {
                    b[j] = r.nextInt(size);
                    a[j] = b[j];
                }
                startTime = System.currentTimeMillis();
                quickSortArbitraryPivot(a);
                endTime = System.currentTimeMillis();
                sumQuick += endTime - startTime;
                startTime = System.currentTimeMillis();
                countingSort(b, size);
                endTime = System.currentTimeMillis();
                sumCounting += endTime - startTime;
            }
            quickTimes[i] = sumQuick / T;
            countingTimes[i] = sumCounting / T;
        }
        Plotter.plot("Counting sort on arrays with elements < n", countingTimes,
                "Quick sort on arrays with elements < n", quickTimes);

    }

    /**
     * Compares the selection sort algorithm against quick sort on random arrays
     */
    public static void insertionVsQuick() {
        double[] quickTimes = new double[INSERTION_VS_QUICK_LENGTH];
        double[] insertionTimes = new double[INSERTION_VS_QUICK_LENGTH];
        long startTime, endTime;
        Random r = new Random();
        for (int i = 0; i < INSERTION_VS_QUICK_LENGTH; i++) {
            long sumQuick = 0;
            long sumInsertion = 0;
            for (int k = 0; k < T; k++) {
                int size = (int) Math.pow(2, i);
                double[] a = new double[size];
                double[] b = new double[size];
                for (int j = 0; j < a.length; j++) {
                    a[j] = r.nextGaussian() * 5000;
                    b[j] = a[j];
                }
                startTime = System.currentTimeMillis();
                quickSortArbitraryPivot(a);
                endTime = System.currentTimeMillis();
                sumQuick += endTime - startTime;
                startTime = System.currentTimeMillis();
                insertionSort(b);
                endTime = System.currentTimeMillis();
                sumInsertion += endTime - startTime;
            }
            quickTimes[i] = sumQuick / T;
            insertionTimes[i] = sumInsertion / T;
        }
        Plotter.plot("quick sort on random array", quickTimes, "insertion sort on random array", insertionTimes);
    }

    /**
     * Compares the merge sort algorithm against quick sort on random arrays
     */
    public static void mergeVsQuick() {
        double[] quickTimes = new double[MERGE_VS_QUICK_LENGTH];
        double[] mergeTimes = new double[MERGE_VS_QUICK_LENGTH];
        long startTime, endTime;
        Random r = new Random();
        for (int i = 0; i < MERGE_VS_QUICK_LENGTH; i++) {
            long sumQuick = 0;
            long sumMerge = 0;
            for (int k = 0; k < T; k++) {
                int size = (int) Math.pow(2, i);
                double[] a = new double[size];
                double[] b = new double[size];
                for (int j = 0; j < a.length; j++) {
                    a[j] = r.nextGaussian() * 5000;
                    b[j] = a[j];
                }
                startTime = System.currentTimeMillis();
                quickSortArbitraryPivot(a);
                endTime = System.currentTimeMillis();
                sumQuick += endTime - startTime;
                startTime = System.currentTimeMillis();
                mergeSort(b);
                endTime = System.currentTimeMillis();
                sumMerge += endTime - startTime;
            }
            quickTimes[i] = sumQuick / T;
            mergeTimes[i] = sumMerge / T;
        }
        Plotter.plot("quick sort on random array", quickTimes, "merge sort on random array", mergeTimes);
    }

    /**
     * Compares the merge sort algorithm against quick sort on pre-sorted arrays
     */
    public static void insertionVsQuickOnSortedArray() {
        double[] quickTimes = new double[INSERTION_VS_QUICK_SORTED_LENGTH];
        double[] insertionTimes = new double[INSERTION_VS_QUICK_SORTED_LENGTH];
        long startTime, endTime;
        for (int i = 0; i < INSERTION_VS_QUICK_SORTED_LENGTH; i++) {
            long sumQuick = 0;
            long sumInsertion = 0;
            for (int k = 0; k < T; k++) {
                int size = (int) Math.pow(2, i);
                double[] a = new double[size];
                double[] b = new double[size];
                for (int j = 0; j < a.length; j++) {
                    a[j] = j;
                    b[j] = j;
                }
                startTime = System.currentTimeMillis();
                quickSortArbitraryPivot(a);
                endTime = System.currentTimeMillis();
                sumQuick += endTime - startTime;
                startTime = System.currentTimeMillis();
                insertionSort(b);
                endTime = System.currentTimeMillis();
                sumInsertion += endTime - startTime;
            }
            quickTimes[i] = sumQuick / T;
            insertionTimes[i] = sumInsertion / T;
        }
        Plotter.plot("quick sort on sorted array", quickTimes, "insertion sort on sorted array", insertionTimes);
    }

    /**
     * Compares the quick sort algorithm once with a choice of an arbitrary pivot
     * and once with a choice of a random pivot
     */
    public static void arbitraryPivotVsRandomPivot() {
        double[] arbitraryTimes = new double[ARBITRARY_VS_RANDOM_LENGTH];
        double[] randomTimes = new double[ARBITRARY_VS_RANDOM_LENGTH];
        long startTime, endTime;
        Random r = new Random();
        for (int i = 0; i < ARBITRARY_VS_RANDOM_LENGTH; i++) {
            long sumArbitrary = 0;
            long sumRandom = 0;
            for (int k = 0; k < T; k++) {
                int size = (int) Math.pow(2, i);
                double[] a = new double[size];
                double[] b = new double[size];
                for (int j = 0; j < a.length; j++) {
                    a[j] = r.nextGaussian() * 5000;
                    b[j] = a[j];
                }
                startTime = System.currentTimeMillis();
                quickSortArbitraryPivot(a);
                endTime = System.currentTimeMillis();
                sumArbitrary += endTime - startTime;
                startTime = System.currentTimeMillis();
                quickSortRandomPivot(b);
                endTime = System.currentTimeMillis();
                sumRandom += endTime - startTime;
            }
            arbitraryTimes[i] = sumArbitrary / T;
            randomTimes[i] = sumRandom / T;
        }
        Plotter.plot("quick sort with an arbitrary pivot", arbitraryTimes, "quick sort with a random pivot",
                randomTimes);
    }

}
