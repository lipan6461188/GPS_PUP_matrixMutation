#!/opt/local/bin/perl -w
use strict;
use 5.010;
use Getopt::Long;       #长参数


#TP : 309
#TN : 3965

# --upper   阈值上界
# --lower   阈值下界
# --gap     阈值gap
# --blosum  BLOSUM矩阵的文件
# --sample  正负样本的文件
# --output  输出文件

#输入序列的长度、阈值下界与上界、gap、blosum矩阵文件、训练集与测试集文件
our $upper = 4;
our $lower = 0;
our $gap = 0.0001;
our $blosum;
our $sample;
our $output;
Getopt::Long::GetOptions('upper=s' => \$upper, 'lower=s' => \$lower, 'gap=s' => \$gap, 'blosum=s' => \$blosum, 'sample=s' => \$sample, 'output=s' => \$output)
or die("Error in command line arguments\n");

if(!defined($blosum) || !defined($sample) || !defined($output)){ warn("请输如sample样本文件和blosum矩阵\n"); exit(-1);}

#权重分数
our @scores = ( 1,1,1,1,1,
                1,1,1,1,1,
                1,1,1,1,1,
                1,1,1,1,1,
                1,1,1
               );
#引入前一步的结果
$scores[4] = 0;
$scores[7] = 0;
$scores[9] = 0;

#读入文件
open(HANDLE, "<$sample") || die "打开文件失败\n";
chomp(my @all_seq = <HANDLE>);
close HANDLE;

#正负样本分类
our @positive = ();
our @negative = ();
our @positive_score = ();
our @negative_score = ();
our $Sn_old = 0;
foreach my $eachSeq(@all_seq)
{
    my @seq_prop = split(/\t/,$eachSeq);
    if($seq_prop[1] eq "+")
    { push(@positive, $seq_prop[0]); }
    else
    { push(@negative, $seq_prop[0]); }
}

say '正样本个数：',$#positive+1;
say '负样本个数: ',$#negative+1;
our $pos_num = $#positive+1;
our $neg_num = $#negative+1;

#获取BLOSUM62矩阵
our %BLOSUM = constructBLOSUM62();
#构造一个偏移矩阵
our %offset = %BLOSUM;

#偏移矩阵所有的值都为零
while(my $key = each %offset)
{
    $offset{$key} = 0;
}


my @performance = &performance();   #首先计算一下分数，Sn_old
say "BEGIN ",pop(@performance),"\t",pop(@performance),"\t",pop(@performance);

&matrixMutation();

&writeToFile();

exit(0);

sub writeToFile
{
    open OUTPUT ">$output" or die("无法打开输出文件\n");
    
}

#矩阵突变
sub matrixMutation
{
    my @keys = keys %BLOSUM;
    my @values = values %BLOSUM;
    my $blosum_size = scalar(@keys);
    
    for(my $index=0; $index < 10000; $index++)
    {
        #随机地取出一个位置做突变
        my $rad = int(rand($blosum_size));
        
        $offset{$keys[$rad]}++;
        my @performance = &performance();
        if( !pop(@performance) )
        {
            say "- $keys[$rad] +1 ",pop(@performance),"\t",pop(@performance);
            $offset{$keys[$rad]} -= 2;
            #一个位点减1以后性能还是下降,就跳过这个位点
            @performance = &performance();
            if(!pop(@performance))
            {
                say "- $keys[$rad] -1 ",pop(@performance),"\t",pop(@performance);
                $offset{$keys[$rad]}++;  #回到原点
            }else{
                say "+ $keys[$rad] -1 ",pop(@performance),"\t",pop(@performance);
            }
        }
        else{
            say "+ $keys[$rad] +1 ",pop(@performance),"\t",pop(@performance);
        }
    }
    
}

#看看矩阵修改以后，性能有没有提高
#返回值是Sp，Sn和是否提高的标志
# @items = {Sp, Sn, flag}
sub performance
{
    my $lower = $lower;
    my $upper = $upper;
    my $gap = $gap;
    
    #首先计算所有样本的分数
    &computeScoreForEveryElem();
    my @items = &iter_compute_Sn($lower,$upper,$gap);
    while( abs($items[0] - 0.91) > 0.01 )
    {
        $gap /= 10;
        @items = &iter_compute_Sn($lower,$upper,$gap);
    }
    pop(@items);
    return @items;
}

#计算Sn和Sp的函数,参数是阈值
sub compute_Sn_Sp
{
    my $threshold = shift @_;
    
    #初始化这些变量
    my $TP=0;
    my $FN=0;
    my $TN=0;
    my $FP=0;
    my $Sn=0;
    my $Sp=0;
    
    foreach(@positive_score)
    {
        if($_ >= $threshold)
        { $TP++; }
        else
        { $FN++; }
    }
    
    foreach(@negative_score)
    {
        if($_ < $threshold)
        { $TN++; }
        else
        { $FP++; }
    }
    
    $Sp = $TP/($TP+$FN);
    $Sn = $TN/($TN+$FP);
    
    return ($Sn,$Sp);
}

#参数分别是阈值下界，上界和gap
sub iter_compute_Sn
{
    my $lower = shift @_;
    my $upper = shift @_;
    my $gap = shift @_;
    
    for(my $threshold=$lower; $threshold <= $upper; $threshold += $gap)
    {
        my @Sn_Sp = &compute_Sn_Sp($threshold);
        
        #定义返回项.返回项中分别是Sp, Sn, 是否提高, 阈值
        my @return_items;
        push(@return_items,$Sn_Sp[1]);
        push(@return_items,$Sn_Sp[0]);
        if( $Sn_Sp[1] - 0.91 <= 0.01 || $Sn_Sp[1] < 0.91)
        {
            if( $Sn_old < $Sn_Sp[0] ) #性能提高了
            {
                $Sn_old = $Sn_Sp[0];
                push(@return_items,1);
            }else{              #性能没有提高
                push(@return_items,0);
            }
            push(@return_items,$threshold);
            return @return_items;
        }
        
    }
}


#计算每一个样本的分数。返回正样本平均分数、负样本平均分数、差值
sub computeScoreForEveryElem
{
    #=============计算正样本==============#
    my $pos_s = 0;
    @positive_score = ();
    foreach my $one (@positive)
    {
        #  say $temp;
        my $sum = 0;#定义一个和，最后要求平均。
        foreach my $each (@positive)
        {
            #不要求自己与自己的相似度，因为这个结果会很高
            if( $one eq  $each) { next; }
            #求两条序列之间的距离
            my $s_grade = similarity($one, $each);
            #如果序列的距离小于0，就设置为0，只把分数大于0的序列分数相加
            if($s_grade > 0) { $sum += $s_grade;}
        }
        #say "Ave: ".$average_grade;
        #求一条序列与其他序列距离的平均值，把这个平均值放置到数组中，今后用于统计
        push(@positive_score, $sum/($#positive+1));
        $pos_s += $sum;
        # say "正样本$one的得分是: ".$sum;
    }
    
    #=============计算负样本==============#
    @negative_score = ();
    my $neg_s = 0;
    foreach my $one (@negative)
    {
        my $sum = 0;#定义一个和，最后要求平均。
        foreach my $each (@positive)
        {
            #求两条序列之间的距离
            my $s_grade = similarity($one, $each);
            #如果序列的距离小于0，就设置为0
            if($s_grade > 0) { $sum += $s_grade;}
        }
        #求一条序列与其他序列距离的平均值，把这个平均值放置到数组中，今后用于统计
        push(@negative_score, $sum/($#positive+1));
        $neg_s += $sum;
        # say "负样本$one的得分是: ".$sum;
    }
    
    my $overlap = $pos_s/$#positive - $neg_s/$#negative;
    return ($pos_s/$#positive, $neg_s/$#negative, $overlap);
}


#子函数用于构造哈希
sub constructBLOSUM62
{
    #读入哈希
    open HASH,"<$blosum";
    chomp(my @file = <HASH>);
    close HASH;
    
    my %BLOSUM62;
    my @index;
    
    foreach my $eachline(@file)
    {
        #截取信息
        my @info = split(/ +/,$eachline);
        #第一行的信息
        if($info[0] eq "")
        {
            shift @info;
            @index = @info;
            next;
        }
        #其余行的信息
        my $tablet = shift @info;
        for(my $i=0; $i < 24; $i++)
        {
            $BLOSUM62{ $tablet.$index[$i] } = $info[$i];
        }
    }
    return %BLOSUM62;
}

#用于确定两条序列的相似度
sub similarity
{
    my $seq1 = pop @_;
    my $seq2 = pop @_;
    if(length($seq1) ne length($seq2))
    {
        say '$seq1: ',$seq1;
        say '$seq2: ',$seq2;
        warn("==>要比较的序列长度不等<==\n");
        exit(0);
    }
    #  say "==============";
    my $score = 0;
    for (my $var = 0; $var < length($seq1); $var++)
    {
        my $bp1 = substr ($seq1,$var,1);
        my $bp2 = substr ($seq2,$var,1);
        $score += ($BLOSUM{ $bp1.$bp2 }+$offset{ $bp1.$bp2 }) * $scores[$var];
        #    say $weights[ $var ];
    }
    
    return $score;
}





