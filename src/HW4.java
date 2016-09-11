import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Stack;

public class HW4 {

  // public static enum Command {
  // L, R, INC, DEC, OB, CB, IN, OUT
  // }

  // static final Command[] commandlist = { Command.L, Command.R, Command.INC,
  // Command.DEC, Command.OB, Command.CB, Command.IN, Command.OUT };

  static final int L = 1;
  static final int R = 2;
  static final int INC = 3;
  static final int DEC = 4;
  static final int OB = 5;
  static final int CB = 6;
  static final int IN = 7;
  static final int OUT = 8;

  public static class CommandPair {
    int type;
    int run;

    public CommandPair(int com, int run) {
      this.type = com;
      this.run = run;
    }

    public String toString() {
      String res = "";
      String typeStr = commandToString(type);
      for (int i = 0; i < run; i++) {
        res += typeStr;
      }
      return res;
    }
  }

  public static String commandToString(int type) {
    switch (type) {
    case L:
      return "<";
    case R:
      return ">";
    case INC:
      return "+";
    case DEC:
      return "-";
    case OB:
      return "[";
    case CB:
      return "]";
    case IN:
      return ",";
    case OUT:
      return ".";
    }
    return "?";
  }

  public static class CommandPairList {
    private ArrayList<CommandPair> coms = new ArrayList<CommandPair>();
    int lengthSum = 0;

    public void add(CommandPair com) {
      coms.add(com);
      lengthSum += com.run;
    }

    public void add(int loc, CommandPair com) {
      coms.add(loc, com);
      lengthSum += com.run;
    }

    public void addAll(CommandPairList data) {
      this.addAll(this.coms.size(), data);
    }

    public void addAll(int loc, CommandPairList data) {
      coms.addAll(loc, data.coms);
      lengthSum += data.lengthSum;
    }

    public void addAll(int[] data) {
      this.addAll(this.coms.size(), data);
    }

    public void addAll(int loc, int[] data) {
      CommandPairList toAdd = new CommandPairList();
      for (int i = 0; i < data.length; i++) {
        if (toAdd.coms.size() == 0
            || toAdd.get(toAdd.coms.size() - 1).type != data[i]) {
          toAdd.add(new CommandPair(data[i], 0));
        }
        toAdd.get(toAdd.coms.size() - 1).run++;
        toAdd.lengthSum++;
      }
      this.addAll(loc, toAdd);
    }

    public CommandPair get(int loc) {
      return coms.get(loc);
    }

    public String toString() {
      String res = lengthSum + ": ";
      for (CommandPair com : coms) {
        res += com;
      }
      return res;
    }

    public int[] toCommandArray() {
      int[] res = new int[lengthSum];
      int i = 0;
      for (CommandPair com : coms) {
        for (int j = 0; j < com.run; j++, i++) {
          res[i] = com.type;
        }
      }
      return res;
    }
  }

  public static final byte[] targetString = stringToByteArray("Hello, world!");

  static PrintWriter out;

  public static void main(String[] args) throws FileNotFoundException {
    out = new PrintWriter(new File("outHW4.txt"));
    // TODO
    // byte[] target = stringToByteArray("Hello, World!");
    // CommandPairList prog = simplePathfinder(new byte[15], 5, target, 0,
    // Integer.MAX_VALUE);
    // System.out.println(eval(prog.toCommandArray()));
    // cycleLengths(15, 30);
    // System.out.println(isValidRecRel(new int[] { DEC, DEC, L, DEC, L, L, INC,
    // OB, INC, OB, L, INC, R, DEC, DEC, DEC, R, DEC, R, DEC, R, DEC, L, L, L,
    // CB, R, CB }));
    // printNestedLoops2("abcdefghijklmnopqr");
    for (int i = 1; i < 50; i++) {
      iterRecRel4(i);
    }

    // accum[precon[<lhmult>deltas<<<]lhadd>postcon]
    System.out.println();
    TapeState benchmark = evalRecRel(new byte[] { 1, 0, -1, -2 }, (byte) 1,
        (byte) 1, new byte[] { -3, -1, -1, -1 }, (byte) 0, (byte) 0);
    System.out.println(benchmark);
    System.out.println(simplePathfinder(benchmark.tapeArray(), benchmark.dp,
        targetString, 0, Integer.MAX_VALUE));
    System.out
        .println(progToString(recRelToProg(new byte[] { 1, 0, -1, -2 },
            (byte) 1, (byte) 1, new byte[] { -3, -1, -1, -1 }, (byte) 0,
            (byte) 0)));
    System.out.println(256 % 3);
  }

  public static byte[] stringToByteArray(String s) {
    byte[] res = new byte[s.length()];
    for (int i = 0; i < s.length(); i++) {
      res[i] = (byte) s.charAt(i);
    }
    return res;
  }

  static CommandPairList simplePathfinder(final byte[] tape, final int loc,
      final byte[] reqout, final int curchar, int best) {
    CommandPairList bestSuffix = null;
    for (Integer i = loc; i != null && i > tape.length - 15; i = staggerNext(i,
        loc, tape.length)) {
      CommandPairList curSuffix = new CommandPairList();
      if (i < loc) {
        curSuffix.add(new CommandPair(L, loc - i));
      } else if (i > loc) {
        curSuffix.add(new CommandPair(R, i - loc));
      }
      int dif = (byte) (reqout[curchar] - tape[i]);
      int absdif = dif;
      if (absdif < 0) {
        absdif = -absdif;
      }
      if (dif > 0) {
        curSuffix.add(new CommandPair(INC, absdif));
      } else if (dif < 0) {
        curSuffix.add(new CommandPair(DEC, absdif));
      }
      curSuffix.add(new CommandPair(OUT, 1));
      if (curSuffix.lengthSum >= best) {
        continue;
      }
      if (curchar < reqout.length - 1) {
        byte[] childtape = tape.clone();
        childtape[i] = reqout[curchar];
        CommandPairList childSuf = simplePathfinder(childtape, i, reqout,
            curchar + 1, best - curSuffix.lengthSum);
        if (childSuf == null) {
          continue;
        }
        curSuffix.addAll(childSuf);
      }
      // if (curchar == 0) {
      // System.out.println(curSuffix);
      // }
      best = curSuffix.lengthSum;
      bestSuffix = curSuffix;
    }
    return bestSuffix;
  }

  public static Integer staggerNext(int cur, int pivot, int limit) {
    Integer res = null;
    if (cur > pivot) {
      res = (pivot << 1) - cur;
      if (res < 0) {
        res = cur + 1;
        if (res >= limit) {
          return null;
        }
      }
    } else {
      res = (pivot << 1) - cur + 1;
      if (res >= limit) {
        res = cur - 1;
        if (res < 0) {
          return null;
        }
      }
    }
    return res;
  }

  public static class TapeState {
    ArrayList<Byte> tape;
    ArrayList<Character> output;
    int dp;
    int time;

    public TapeState() {
      tape = new ArrayList<Byte>();
      output = new ArrayList<Character>();
      dp = 0;
      time = 0;
    }

    public String toString() {
      return output + "" + tape + dp;
    }

    public byte[] tapeArray() {
      byte[] res = new byte[tape.size()];
      for (int i = 0; i < tape.size(); i++) {
        res[i] = tape.get(i);
      }
      return res;
    }

    public byte[] toByteArray() {
      // byte[] res = new byte[8+output.size()+tape.size()];
      byte[] res = new byte[0 + tape.size()];
      int j = 0;
      // res[j++] = (byte) ((dp >>> 0) & 0xff);
      // res[j++] = (byte) ((dp >>> 8) & 0xff);
      // res[j++] = (byte) ((dp >>> 16) & 0xff);
      // res[j++] = (byte) ((dp >>> 24) & 0xff);
      // res[j++] = (byte) ((output.size() >>> 0) & 0xff);
      // res[j++] = (byte) ((output.size() >>> 8) & 0xff);
      // res[j++] = (byte) ((output.size() >>> 16) & 0xff);
      // res[j++] = (byte) ((output.size() >>> 24) & 0xff);
      // for(int i = 0; i < output.size(); i++){
      // res[j++] = (byte)(output.get(i).charValue());
      // }
      int start = 0;
      if (tape.size() > 15) {
        start = tape.size() - 15;
      }
      for (int i = start; i < tape.size(); i++) {
        res[j++] = tape.get(i);
      }
      return res;
    }
  }

  public static void printNestedLoops2(String s) {
    for (int i = 0; i < s.length(); i++) {
      char curchar = s.charAt(i);
      if (i == 0) {
        System.out.println("for(int " + curchar + "x = -maxsum; " + curchar
            + "x <= maxsum; " + curchar + "x++){");
        System.out.println("int sum" + curchar + "x = maxsum - Math.abs("
            + curchar + "x);");
      } else {
        char prevchar = s.charAt(i - 1);
        System.out.println("for(int " + curchar + "x = -sum" + prevchar + "x; "
            + curchar + "x <= sum" + prevchar + "x; " + curchar + "x++){");
        System.out.println("int sum" + curchar + "x = sum" + prevchar
            + "x - Math.abs(" + curchar + "x);");
      }
    }
    for (int i = 0; i < s.length(); i++) {
      System.out.print("}");
    }
    System.out.println();
  }

  public static void iterRecRel4(int maxsum) {

    int bfsize = 50000000;
    BloomFilter<TapeState> bf = new BloomFilter<TapeState>(0.01, bfsize);
    long dupecounter = 0;
    long bfcounter = 0;
    long infloopcounter = 0;
    long suboptcounter = 0;
    CommandPairList bestResult = null;
    int bestSuffixLength = Integer.MAX_VALUE;
    boolean bfwarning = false;

    for (int ax = 1; ax <= maxsum; ax++) {
      int sumax = maxsum - Math.abs(ax);
      for (int bx = -sumax; bx <= sumax; bx++) {
        int sumbx = sumax - Math.abs(bx);
        // for (int cx = -sumbx; cx <= sumbx; cx++) {
        // int sumcx = sumbx - Math.abs(cx);
        // for (int dx = -sumcx; dx <= sumcx; dx++) {
        // int sumdx = sumcx - Math.abs(dx);
        // for (int ex = -sumdx; ex <= sumdx; ex++) {
        // int sumex = sumdx - Math.abs(ex);
        // for (int fx = -sumex; fx <= sumex; fx++) {
        // int sumfx = sumex - Math.abs(fx);
        // for (int gx = -sumfx; gx <= sumfx; gx++) {
        // int sumgx = sumfx - Math.abs(gx);

        byte[] accum = { (byte) ax, (byte) bx, (byte) 0, (byte) 0, (byte) 0,
            (byte) 0, (byte) 0 };
        for (int hx = -sumbx; hx <= sumbx; hx++) {
          int sumhx = sumbx - Math.abs(hx);

          byte precon = (byte) hx;
          for (int ix = -sumhx; ix <= sumhx; ix++) {
            int sumix = sumhx - Math.abs(ix);

            if (ix != 0) {
              byte lhmult = (byte) ix;
              for (int jx = 1; jx <= 7; jx += 2) {
                int sumjx = sumix - Math.abs(jx);
                for (int kx = -sumjx; kx <= sumjx; kx++) {
                  int sumkx = sumjx - Math.abs(kx);
                  for (int lx = -sumkx; lx <= sumkx; lx++) {
                    int sumlx = sumkx - Math.abs(lx);
                    for (int mx = -sumlx; mx <= sumlx; mx++) {
                      int summx = sumlx - Math.abs(mx);
                      for (int nx = -summx; nx <= summx; nx++) {
                        int sumnx = summx - Math.abs(nx);
                        for (int ox = -sumnx; ox <= sumnx; ox++) {
                          int sumox = sumnx - Math.abs(ox);
                          for (int px = -sumox; px <= sumox; px++) {
                            // if (px != 0 || ox != 0) {
                            int sumpx = sumox - Math.abs(px);

                            byte[] deltas = { (byte) jx, (byte) kx, (byte) lx,
                                (byte) mx, (byte) nx, (byte) ox, (byte) px };
                            for (int qx = 0; qx <= 0; qx++) {
                              int sumqx = sumpx - Math.abs(qx);
                              byte lhadd = (byte) qx;
                              for (int rx = -sumqx, flipflop = 0; rx <= sumqx
                                  && flipflop < 2; rx += sumqx, flipflop++) {
                                if (rx == 0 && flipflop == 1) {
                                  break;
                                }
                                int sumrx = sumqx - Math.abs(rx);
                                byte postcon = (byte) rx;

                                // System.out.println(progToString(recRelToProg(
                                // accum.clone(), precon, lhmult,
                                // deltas.clone(), lhadd,
                                // postcon)));

                                TapeState res = evalRecRel(accum.clone(),
                                    precon, lhmult, deltas, lhadd, postcon);

                                if (res == null) {
                                  infloopcounter++;
                                } else {
                                  if (bf.contains(res.toByteArray())) {
                                    dupecounter++;
                                  } else {
                                    bfcounter++;
                                    bf.add(res.toByteArray());
                                    if (!bfwarning && bfcounter > bfsize) {
                                      bfwarning = true;
                                      out.println("=========BLOOM FILTER OVERFLOW==========");
                                      out.flush();
                                    }

                                    CommandPairList suffix = simplePathfinder(
                                        res.tapeArray(), res.dp, targetString,
                                        0, bestSuffixLength - maxsum + sumrx
                                            - 18 + 1);
                                    if (suffix == null) {
                                      suboptcounter++;
                                    } else {
                                      bestResult = suffix;
                                      int[] prog = recRelToProg(accum.clone(),
                                          precon, lhmult, deltas.clone(),
                                          lhadd, postcon);
                                      bestResult.addAll(0, prog);
                                      // if (bestResult.lengthSum < 85) {
                                      out.println(res);
                                      out.println(bestResult);
                                      out.println("len:" + (maxsum - sumrx)
                                          + " unique:" + bfcounter + " dupe:"
                                          + dupecounter + " inf:"
                                          + infloopcounter + " subopt:"
                                          + suboptcounter);
                                      out.flush();
                                      // }
                                      bestSuffixLength = bestResult.lengthSum;
                                    }
                                  }
                                }

                              }
                            }
                            // }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
        // }
        // }
        // }
        // }
        // }
      }
    }
  }

  public static int[] recRelToProg(byte[] accum, byte precon, byte lhmult,
      byte[] deltas, byte lhadd, byte postcon) {
    ArrayList<Integer> prog = new ArrayList<Integer>();
    int proglength = 25 + Math.abs(precon) + Math.abs(lhmult) + Math.abs(lhadd)
        + Math.abs(postcon);
    for (int i = 0; i < accum.length; i++) {
      proglength += Math.abs(accum[i]);
    }
    for (int i = 0; i < deltas.length; i++) {
      proglength += Math.abs(deltas[i]);
    }

    boolean posacc = false;
    for (int accpt = accum.length - 1; accpt >= 0; accpt--) {
      if (posacc) {
        prog.add(L);
      } else if (accum[accpt] == 0) {
        continue;
      }
      posacc = true;
      for (int ind = 0; ind < Math.abs(accum[accpt]); ind++) {
        if (accum[accpt] > 0) {
          prog.add(INC);
        } else {
          prog.add(DEC);
        }
      }
    }

    prog.add(OB);

    for (int ind = 0; ind < Math.abs(precon); ind++) {
      if (precon > 0) {
        prog.add(INC);
      } else {
        prog.add(DEC);
      }
    }

    prog.add(OB);
    prog.add(L);

    for (int ind = 0; ind < Math.abs(lhmult); ind++) {
      if (lhmult > 0) {
        prog.add(INC);
      } else {
        prog.add(DEC);
      }
    }

    int lastnonzerodelta = 0;
    for (int i = 1; i < deltas.length; i++) {
      if (deltas[i] != 0) {
        lastnonzerodelta = i;
      }
    }

    for (int delpt = 0; delpt <= lastnonzerodelta; delpt++) {
      prog.add(R);
      for (int ind = 0; ind < Math.abs(deltas[delpt]); ind++) {
        if (deltas[delpt] > 0) {
          prog.add(INC);
        } else {
          prog.add(DEC);
        }
      }
    }

    for (int delpt = 0; delpt < lastnonzerodelta; delpt++) {
      prog.add(L);
    }

    prog.add(CB);

    for (int ind = 0; ind < Math.abs(lhadd); ind++) {
      if (lhadd > 0) {
        prog.add(INC);
      } else {
        prog.add(DEC);
      }
    }

    prog.add(R);

    for (int ind = 0; ind < Math.abs(postcon); ind++) {
      if (postcon > 0) {
        prog.add(INC);
      } else {
        prog.add(DEC);
      }
    }

    prog.add(CB);

    int[] res = new int[prog.size()];
    for (int i = 0; i < res.length; i++) {
      res[i] = prog.get(i);
    }

    return res;
  }

  public static String progToString(int[] prog) {
    String res = "";
    res += "{";
    for (int ip = 0; ip < prog.length; ip++) {
      res += commandToString(prog[ip]);
    }
    res += "}";
    return res;
  }

  public static TapeState evalRecRel(byte[] accum, byte precon, byte lhmult,
      byte[] deltas, byte lhadd, byte postcon) {
    // accum[precon[<lhmult>deltas<<<]lhadd>postcon]
    int absd0 = Math.abs(deltas[0]);

    int modinverse = 1;
    switch (absd0) {
    case 3:
      modinverse = 171;
      break;
    case 5:
      modinverse = 205;
      break;
    case 7:
      modinverse = 183;
      break;
    default:
      return null;
    }
    if (deltas[0] < 0) {
      modinverse = -modinverse;
    }

    if (deltas.length != accum.length) {
      return null;
    }
    TapeState tapestate = new TapeState();
    ArrayList<Byte> tape = tapestate.tape;
    // BasicBloomFilter visitedStates = new BasicBloomFilter(0.01, 1000);
    byte lefthand = 0;
    while (tape.size() < 400 && accum[0] != 0) {
      // visitedStates.add(accum);
      accum[0] += precon;
      int factor = -accum[0] * modinverse;
      for (int i = 0; i < deltas.length; i++) {
        accum[i] = (byte) (accum[i] + factor * deltas[i]);
      }
      lefthand = (byte) (lefthand + factor * lhmult);

      tape.add(lefthand);
      lefthand = lhadd;

      for (int i = 0; i < accum.length - 1; i++) {
        accum[i] = accum[i + 1];
      }
      accum[accum.length - 1] = 0;

      accum[0] += postcon;

    }

    if (accum[0] != 0) {
      return null;
    }

    tape.add(lefthand);
    for (int i = 0; i < accum.length - 1; i++) {
      tape.add(accum[i]);
    }

    tapestate.dp = tape.size() - accum.length + 1;

    return tapestate;
  }

}
