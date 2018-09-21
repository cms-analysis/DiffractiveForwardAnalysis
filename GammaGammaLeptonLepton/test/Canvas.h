#ifndef Canvas_h
#define Canvas_h

#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TObjArray.h"
#include "TObjString.h"

#include <string.h>

class PaveText : public TPaveText
{
  public:
    inline PaveText( float x1, float y1, float x2 = kInvalid, float y2 = kInvalid ) :
      TPaveText( x1, y1, ( x2 != kInvalid ) ? x2 : x1+0.25, ( y2 != kInvalid ) ? y2 : y1+0.15, "NB NDC" ) {
      TPaveText::SetFillStyle( 0 );
      TPaveText::SetFillColor( 0 );
      TPaveText::SetLineWidth( 0 );
      TPaveText::SetLineStyle( 0 );
      TPaveText::SetTextFont( 42 );
      TPaveText::SetTextSize( 0.033 );
      TPaveText::SetTextAlign( kHAlignRight+kVAlignBottom );
    }

  private:
    static constexpr float kInvalid = -999.;
};

class Canvas : public TCanvas
{
  public:
    inline Canvas( const char* name, const char* title = "", bool ratio = false ) :
      TCanvas( name, "", 600, 600 ), title_( title ),
      legend_x_pos_( 0.5 ), legend_y_pos_( 0.75 ),
      legend_x_size_( 0.35 ), legend_y_size_( 0.15 ),
      ratio_( ratio ), divided_( false )
    {
      Build();
    }

    inline ~Canvas() {}

    inline void Prettify( TH1* obj ) {
      TAxis* x = dynamic_cast<TAxis*>( obj->GetXaxis() ),
            *y = dynamic_cast<TAxis*>( obj->GetYaxis() );
      x->SetLabelFont( 43 );
      x->SetTitleFont( 43 );
      y->SetLabelFont( 43 );
      y->SetTitleFont( 43 );
      if ( !divided_ ) {
        x->SetLabelSize( 20 );
        x->SetTitleSize( 26 );
        y->SetLabelSize( 20 );
        y->SetTitleSize( 26 );
        y->SetTitleOffset( 1.4 );
      }
      else {
        x->SetLabelSize( 16 );
        x->SetTitleSize( 20 );
        y->SetLabelSize( 15 );
        y->SetTitleSize( 20 );
        x->SetTitleOffset( 1.8 );
        x->SetLabelOffset( 0.02 );
        y->SetTitleOffset( 2.1 );
      }
      x->SetTitleColor( kBlack );
      if ( ratio_ ) {
        x->SetTitleOffset( 3.2 );
        x->SetLabelOffset( 0.02 );
      }

      // axis titles
      TString ttle = obj->GetTitle();
      if ( ttle.Contains( "\\" ) ) {
        TObjArray* tok = ttle.Tokenize( "\\" );
        TString x_title = "", y_title = "", unit = "", form_spec = "";
        if ( tok->GetEntries() > 0 )
          x_title = dynamic_cast<TObjString*>( tok->At( 0 ) )->String();
        if ( tok->GetEntries() > 1 )
          y_title = dynamic_cast<TObjString*>( tok->At( 1 ) )->String();
        if ( tok->GetEntries() > 2 ) {
          unit = dynamic_cast<TObjString*>( tok->At( 2 ) )->String();
          if ( unit.Contains( "?" ) ) { // extract format specifier
            TObjArray* tok2 = unit.Tokenize( "?" );
            if ( tok2->GetEntries()>1 ) {
              unit = dynamic_cast<TObjString*>( tok2->At( 0 ) )->String();
              form_spec = dynamic_cast<TObjString*>( tok2->At( 1 ) )->String();
            }
            else {
              unit = "";
              form_spec = dynamic_cast<TObjString*>( tok2->At( 0 ) )->String();
            }
          }
          if ( !unit.IsNull() or !form_spec.IsNull() ) {
            if ( !unit.IsNull() ) x_title = Form( "%s (%s)", x_title.Data(), unit.Data() );
            if ( !form_spec.IsNull() ) {
              TString format = Form( "%%s / %%%s %%s", form_spec.Data() );
              y_title = Form( format.Data(), y_title.Data(), GetBinning( obj ), unit.Data() );
            }
            else y_title = Form( "%s / %g %s", y_title.Data(), GetBinning( obj ), unit.Data() );
          }
        }
        obj->GetXaxis()->SetTitle( x_title );
        obj->GetYaxis()->SetTitle( y_title );
        obj->SetTitle( "" );
      }
    }

    inline void Divide( int num_cols, int num_lines = 1, float xmargin = 0.01, float ymargin = 0.01, int color = 0 ) override {
      if ( ratio_ )
        throw std::runtime_error( "Canvas is already divided!" );
      TCanvas::Divide( num_cols, num_lines, xmargin, ymargin, color );
      double top_margin = 0.055;
      double pad_x = 1./num_cols, pad_y = ( 1.-top_margin )/num_lines;
      for ( unsigned short l = 0; l < num_lines; ++l ) {
        for ( unsigned short c = 0; c < num_cols; ++c ) { // fetch one line by one line
          TPad* p = dynamic_cast<TPad*>( TCanvas::GetPad( 1+l*num_cols+c ) );
          p->SetPad( 0.+pad_x*c, 1.-top_margin-pad_y*( l+1 ), 0.+pad_x*( c+1 ), 1.-top_margin-pad_y*l );
          p->SetLeftMargin( 0.15 );
          if ( c == num_cols-1 )
            p->SetRightMargin( TCanvas::GetRightMargin() );
          p->SetTopMargin( 0. );
          p->SetBottomMargin( 0.15 );
          p->SetTicks( 1, 1 );
        }
      }
      divided_ = true;
    }

    typedef std::vector< std::pair<std::string,TH1*> > HistsMap;
    inline void RatioPlot( HistsMap hm, float ymin = kInvalid, float ymax = kInvalid, float xline = kInvalid, bool pull = false ) {
      if ( !ratio_ )
        throw std::runtime_error( "Trying to produce a ratio plot with an unsplitted Canvas!" );
      TH1* denom = hm.begin()->second, *numer = nullptr;
      denom->GetXaxis()->SetTitle( "" );
      TCanvas::cd( 2 );
      const char* ytitle = ( pull ) ? "Pull" : "Ratio";

      TH1D* denom_err = (TH1D*)denom->Clone(), *denom_err2 = (TH1D*)denom->Clone();
      denom_err2->Sumw2( false );
      denom_err->Divide( denom_err2 );

      unsigned short i = 0;
      for ( HistsMap::const_iterator it = hm.begin()+1; it != hm.end(); ++it ) {
        numer = dynamic_cast<TH1*>( it->second->Clone() );
        //ratio1->Sumw2(); ratio2->Sumw2();
        if ( pull && !numer->Add( denom, -1. ) )
          throw std::runtime_error( "Failed to add the two histograms before computing pull!" );
        numer->Divide( denom );
        numer->SetLineColor( kBlack );
        numer->SetLineStyle( 1 );
        numer->Draw( ( i++ == 0 ) ? "e" : "e same" );
        //numer->Draw( "p same" );
        if ( ymin != ymax )
          numer->GetYaxis()->SetRangeUser( ymin, ymax );
        Prettify( numer );
        numer->GetYaxis()->SetTitle( Form( "%s%s", ytitle, ( hm.size() > 2 ) ? "s" : "" ) );
      }
      denom_err->Draw( "e2same" );
      denom_err->SetFillColor( kBlack );
      denom_err->SetFillStyle( 3004 );

      if ( xline != kInvalid ) {
        TLine* l = new TLine( denom->GetXaxis()->GetXmin(), xline, denom->GetXaxis()->GetXmax(), xline );
        l->SetLineColor( kBlack );
        l->SetLineWidth( 1 );
        l->Draw();
      }
      Prettify( denom );
      TCanvas::cd();
    }

    inline void SetTopLabel( const char* lab = "" ) {
      TCanvas::cd();
      if ( strcmp( lab, "" ) != 0 )
        title_ = lab;
      if ( !top_label_ )
        BuildTopLabel();
      else
        top_label_->Clear();
      top_label_->AddText( title_ );
    }

    inline TLegend* GetLegend() { return legend_.get(); }
    inline void SetLegendX1( double x ) { legend_x_pos_ = x; }
    inline void SetLegendY1( double y ) { legend_y_pos_ = y; }
    inline void SetLegendSizeX( double x ) { legend_x_size_ = x; }
    inline void SetLegendSizeY( double y ) { legend_y_size_ = y; }
    inline void AddLegendEntry( const TObject* obj, const char* title, Option_t* option = "lpf" ) {
      if ( !legend_ )
        CreateLegend();
      legend_->AddEntry( obj, title, option );
      const unsigned int num_entries = legend_->GetNRows();
      if ( num_entries > 3 ) {
        legend_y_size_ += ( num_entries-3 )*0.015;
        SetLegendY1( legend_->GetY1NDC() );
      }
    }

    inline void Save( const char* ext, const char* out_dir = "." ) {
      TCanvas::cd();
      if ( legend_ && TCanvas::FindObject( legend_.get() ) == 0 )
        legend_->Draw();
      if ( top_label_ && TCanvas::FindObject( top_label_.get() ) == 0 )
        top_label_->Draw();

      const TString ext_str( ext );
      TObjArray* tok = TString( ext ).Tokenize( "," );
      for ( unsigned short i = 0; i < tok->GetEntries(); ++i ) {
        const TString ext_str = dynamic_cast<TObjString*>( tok->At( i ) )->String();
        if ( ext_str == "pdf"
          || ext_str == "png"
          || ext_str == "eps"
          || ext_str == "C" )
          TCanvas::SaveAs( Form( "%s/%s.%s", out_dir, TCanvas::GetName(), ext_str.Data() ) );
      }
    }

  private:
   inline void Build() {
      TCanvas::SetLeftMargin( 0.14 );
      TCanvas::SetTopMargin( 0.06 );
      TCanvas::SetRightMargin( 0.1 );
      TCanvas::SetBottomMargin( 0.12 );
      TCanvas::SetTicks( 1, 1 );

      SetTopLabel();
      if ( ratio_ )
        DivideCanvas();
    }

    inline void DivideCanvas() {
      TCanvas::Divide( 1, 2 );
      TPad* p1 = dynamic_cast<TPad*>( TCanvas::GetPad( 1 ) ),
           *p2 = dynamic_cast<TPad*>( TCanvas::GetPad( 2 ) );
      p1->SetPad( 0., 0.3, 1., 1. );
      p2->SetPad( 0., 0.0, 1., 0.3 );
      p1->SetLeftMargin( TCanvas::GetLeftMargin() );
      p1->SetRightMargin( TCanvas::GetRightMargin() );
      p2->SetLeftMargin( TCanvas::GetLeftMargin() );
      p2->SetRightMargin( TCanvas::GetRightMargin() );
      p1->SetTopMargin( TCanvas::GetTopMargin()+0.025 );
      p1->SetBottomMargin( 0.02 );
      p2->SetTopMargin( 0.02 );
      p2->SetBottomMargin( TCanvas::GetBottomMargin()+0.225 );
      p1->SetTicks( 1, 1 );
      p2->SetTicks( 1, 1 );
      p2->SetGrid( 0, 1 );
      TCanvas::cd( 1 );
    }

    inline void BuildTopLabel() {
      TCanvas::cd();
      top_label_.reset( new PaveText( 0.5, 0.95, 0.925, 0.96 ) );
    }

    inline void CreateLegend() {
      if ( legend_ )
        return;
      if ( ratio_ )
        TCanvas::cd();
      legend_.reset( new TLegend( legend_x_pos_, legend_y_pos_, legend_x_pos_+legend_x_size_, legend_y_pos_+legend_y_size_ ) );
      legend_->SetFillStyle( 0 );
      //legend_->SetLineColor( kWhite );
      //legend_->SetLineColor( kGray );
      //legend_->SetLineWidth( 1 );
      legend_->SetLineWidth( 0 );
      legend_->SetTextSize( 0.035 );
      legend_->Draw();
      if ( ratio_ )
        TCanvas::cd( 1 );
    }
    inline float GetBinning( const TH1* h ) const {
      return ( h->GetXaxis()->GetXmax() - h->GetXaxis()->GetXmin() ) / h->GetXaxis()->GetNbins();
    }

    static constexpr float kInvalid = -999.;
    TString title_;
    std::unique_ptr<PaveText> top_label_;
    std::unique_ptr<TLegend> legend_;
    double legend_x_pos_, legend_y_pos_;
    double legend_x_size_, legend_y_size_;
    bool ratio_;
    bool divided_;
};

template<class T=TH1D> T WithOverflow( T* h ) {
  //function to paint the histogram h with an extra bin for overflows
  unsigned short nx = h->GetNbinsX()+1;
  std::vector<double> xbins( nx+1 );
  for ( unsigned short i = 0; i < nx; ++i )
    xbins[i] = h->GetBinLowEdge( i+1 );
  xbins[nx] = xbins[nx-1] + h->GetBinWidth( nx );
  //book a temporary histogram having extra bins for overflows
  T htmp( Form( "%s_withoverflow", h->GetName() ), h->GetTitle(), nx, &xbins[0] );
  htmp.Sumw2();
  //fill the new histogram including the overflows
  for ( unsigned short i = 1; i <= nx; ++i ) {
    htmp.SetBinContent( htmp.FindBin( htmp.GetBinCenter( i ) ), h->GetBinContent( i ) );
    htmp.SetBinError( htmp.FindBin( htmp.GetBinCenter( i ) ), h->GetBinError( i ) );
  }
  htmp.SetBinContent( htmp.FindBin( h->GetBinLowEdge( 1 )-1 ), h->GetBinContent( 0 ) );
  htmp.SetBinError( htmp.FindBin( h->GetBinLowEdge( 1 )-1 ), h->GetBinError( 0 ) );
  // Restore the number of entries
  htmp.SetEntries( h->GetEffectiveEntries() );
  return htmp;
}

#endif
