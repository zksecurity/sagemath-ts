// Independent oracle: java.util.Random, not a reimplementation of its LCG.
import java.util.Random;

public class RandomReference {
    public static void main(String[] args) {
        Random generator = new Random(Long.parseLong(args[0]));
        int count = Integer.parseInt(args[1]);
        if (args.length > 2 && args[2].equals("bytes")) {
            int skip = Integer.parseInt(args[3]);
            for (int i = 0; i < skip; ++i) generator.nextInt();
            byte[] bytes = new byte[count];
            generator.nextBytes(bytes);
            for (byte value : bytes) System.out.printf("%02x", value & 255);
            System.out.println();
            return;
        }
        for (int i = 0; i < count; ++i) {
            System.out.println(Integer.toUnsignedLong(generator.nextInt()));
        }
    }
}
